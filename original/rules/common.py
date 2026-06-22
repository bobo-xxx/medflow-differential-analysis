from __future__ import annotations
from pathlib import Path
import re
import networkx as nx
import random
import json
import yaml
import shutil
import numpy as np
from .rule import MyRule, RuleSet, Parameter


class VennFiles(MyRule):
    def __init__(self, in_files: dict, venn_plot, out_res, color_panel, confirm_file, prefix):
        super().__init__()
        if isinstance(color_panel, list):
            color_panel = ",".join(color_panel)
        self._input = in_files
        self._output = {"out_res": out_res}
        self._params = {"in_files": Parameter(in_files, value_type="json"), "venn_plot": self.out_prefix_call(venn_plot)}
        self.script = "venn_files.R"
        self.args = ["{params.in_files}", "{params.venn_plot}", "{output.out_res}", color_panel, "{params.confirm_file}", "{params.prefix}"]
        self.out_res = self.out_prefix_call(out_res)

class Report(MyRule):
    def __init__(self, global_confirm, local_confirm, out_word, analysis_name):
        super().__init__()
        self._input = {"global_confirm": global_confirm, "local_confirm": local_confirm}
        self._output = {"out_word": out_word}
        self._params = {"analysis_name": analysis_name}
        self.script = "report.py"
        self.args = ["{input.global_confirm}", "{input.local_confirm}", "{output.out_word}", analysis_name]

def json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return str(obj)
        
# 定义属性转换函数：复杂类型→字符串，基础类型不变
def convert_complex_attrs(value):
    if isinstance(value, (list, dict, set, tuple)):
        # 用 json.dumps() 序列化（比 str() 更规范，支持嵌套结构）
        return json.dumps(value, ensure_ascii=False)
    elif value is None:
        return json.dumps(None)
    else:
        # 基础类型（str、int、float、bool）直接返回
        return value

def revert_complex_attrs(value, original_type=None):
    """
    value: 读取到的属性值（可能是 JSON 字符串）
    original_type: 可选，指定原始类型（如 set/tuple，默认自动推断为 list/dict）
    """
    if isinstance(value, str):  # 仅对字符串尝试反序列化
        try:
            # 先解析为 JSON 原生类型（list/dict）
            parsed = json.loads(value)
            # 根据 original_type 转换为原始类型（如 set/tuple）
            if original_type == set:
                return set(parsed)
            elif original_type == tuple:
                return tuple(parsed)
            else:
                return parsed  # 默认返回 list/dict
        except (json.JSONDecodeError, TypeError):
            # 不是 JSON 字符串，直接返回原始值
            return value
    else:
        # 非字符串类型（int/float/bool），直接返回
        return value

def parse_node(config: dict) -> nx.DiGraph:
    dag = nx.DiGraph()
    
    def check_add_edge(from_node, to_node):
        if from_node not in config:
            raise ValueError(f"{from_node} not found to {to_node}")
        if to_node not in config:
            raise ValueError(f"{to_node} not found")
        dag.add_edge(from_node, to_node)

    for key, value in config.items():
        if key == "step" or key == "global":
            continue
        dag.add_node(key, **value)
        pre = value.get("pre")
        if pre:
            if isinstance(pre, list):
                for p in pre:
                    check_add_edge(p, key)
            elif isinstance(pre, dict):
                for p in pre.values():
                    if isinstance(p, dict):
                        check_add_edge(p["name"], key)
                    elif isinstance(p, str):
                        check_add_edge(p, key)
            else:
                check_add_edge(pre, key)
    if not nx.is_directed_acyclic_graph(dag):
        raise ValueError(f"流程中存在环依赖:\n{list(nx.simple_cycles(dag))}\n无法进行拓扑排序，请检查节点配置！")
    for v, name in enumerate(nx.topological_sort(dag)):
        dag.nodes[name]["index"] = v + 1
    return dag

def compare_dag_graph(old_dag: nx.DiGraph, new_dag: nx.DiGraph, print_diff=True) -> bool:
    """compare_dag_graph.
    Compare two DAGs and return True if they are equal, False otherwise.

    Args:
        old_dag (nx.DiGraph): old_dag
        new_dag (nx.DiGraph): new_dag
        print_diff (bool, optional): if print diff. Defaults to True.

    Returns:
        bool: return True if they are equal, False otherwise.
    """
    old_nodes = set(old_dag.nodes)
    new_nodes = set(new_dag.nodes)
    old_edges = set(old_dag.edges)
    new_edges = set(new_dag.edges)
    added_nodes = new_nodes - old_nodes
    removed_nodes = old_nodes - new_nodes
    added_edges = new_edges - old_edges
    removed_edges = old_edges - new_edges
    if print_diff:
        if added_nodes:
            print("Added nodes:", ", ".join(added_nodes))
        if removed_nodes:
            print("Removed nodes:", ", ".join(removed_nodes))
        if added_edges:
            print("Added edges:", ", ".join(map(str, added_edges)))
        if removed_edges:
            print("Removed edges:", ", ".join(map(str, removed_edges)))
    return len(added_nodes) + len(removed_nodes) + len(added_edges) + len(removed_edges) == 0

def load_dag(dag_path: str | Path) -> nx.DiGraph | None:
    if Path(dag_path).exists():
        dag = nx.read_graphml(dag_path)
        for node, attrs in dag.nodes(data=True):
            for key, value in attrs.items():
                attrs[key] = revert_complex_attrs(value)
        return dag
    else:
        return None

def save_dag(dag: nx.DiGraph, path: str | Path="dag.graphml") -> None:
    dag = dag.copy()
    # 转换节点属性
    for node, attrs in dag.nodes(data=True):
        for key, value in attrs.items():
            attrs[key] = convert_complex_attrs(value)

    # 转换边属性
    for u, v, attrs in dag.edges(data=True):
        for key, value in attrs.items():
            attrs[key] = convert_complex_attrs(value)
    nx.write_graphml(dag, path)

def get_nodes_by_analysis(config: dict, analysis: str) -> set:
    node_set = set()
    for name in config:
        if isinstance(config[name], dict) and config[name].get("analysis", name) == analysis:
            node_set.add(name)
    return node_set

  
class ColorAssigner:
    def __init__(self):
        self.group_colors = [("#84B1ED", "#f199bc"), ("#D09E88", "#dda0dd"),("#FFBE7A","#DE6449"),("#30A9DE","#F68657"),
                             ("#96C6EE","#EEB69A"),("#BCB8D3","#eb998b"),("#DADDFC","#FCBAD3"),("#ffda8e","#ee6e9f"),("#c8c8a9","#F16B6F")]
        self.color_panel = ['#58C9B9','#FFE4E1','#F6B352','#F0F3BD','#eb9f9f','#6a60a9','#FFCF96','#5CAB7D','#754F44','#FDFD96','#9A8C98']
        self.low_colors = ['#4DBBD5','#6AAFE6','#03a6ff']
        self.high_colors = ['#E64B35','#f9320c','#ed317f']
    
    def get_group_color(self):
        colors = random.sample(self.group_colors, 1)
        self.group_colors.remove(colors[0])
        return colors[0]
    
    def get_colors(self, k: int = -1):
        if k == -1 or k > len(self.color_panel):
            k = len(self.color_panel)
        colors = random.sample(self.color_panel, k)
        return colors
    
    def get_low_high_color(self):
        low_colors = random.sample(self.low_colors, 1)
        self.low_colors.remove(low_colors[0])
        high_colors = random.sample(self.high_colors, 1)
        self.high_colors.remove(high_colors[0])
        return low_colors[0], high_colors[0]
    
    def get_heatmap_color(self):
        low_colors = random.sample(self.low_colors, 1)
        high_colors = random.sample(self.high_colors, 1)
        return low_colors[0], "white", high_colors[0]
    
class Cache:
    def __init__(self, config: dict, cahe_path = "cahce", result_path: str | Path = "result"):
        self.cahe_path = cahe_path
        self.result_path = result_path
        self.dag_path = Path(cahe_path, "dag.graphml")
        self.config_bak = Path(cahe_path, "config.yaml")
        self.color_assigner = ColorAssigner()
        self.config = config
        self.load()
        self.assign_color()
        self.save()

    def load(self) -> None:
        if not Path(self.cahe_path).exists():
            Path(self.cahe_path).mkdir(parents=True)
        self.dag = parse_node(config=self.config)
        self.old_dag = load_dag(self.dag_path)
        if Path(self.config_bak).exists():
            self.old_config = yaml.load(open(self.config_bak), Loader=yaml.FullLoader)
            # input_node = self.old_config["input"]
            # input_change = False
            # if input_node["disease_name"] and input_node["gse"]:
            #     if self.config["input"]["disease_name"] != input_node["disease_name"] or self.config["input"]["gse"] != input_node["gse"]:
            #         input_change = True
            # else:
            #     if self.config["input"]["group"] != input_node["group"] or self.config["input"]["data"] != input_node["data"]:
            #         input_change = True
            # if input_change:
            #     print(f"input 模块的输入数据发生改变，是否删除缓存参数（配色等）和结果（{self.result_path}）？")
            common_change = False
            if self.old_dag is not None:
                common_change =  not compare_dag_graph(self.old_dag, self.dag, print_diff=True)
                if common_change:
                    print(f"缓存流程与配置不一致，是否删除缓存参数（配色等）和结果（{self.result_path}）？")
            if common_change:
                s = input("输入 y 删除缓存参数，其他任意输入默认保留：")
                if s == "y":
                    self.old_config = None
                    self.old_dag = None
                    if Path(self.result_path).exists():
                        shutil.rmtree(self.result_path)
                    # raise ValueError("缓存流程与配置不兼容，请检查配置文件或者删除缓存文件：{}".format(self.dag_path))
        else:
            self.old_config = None

    def save(self) -> None:
        save_dag(self.dag, self.dag_path)
        yaml.dump(self.config, open(self.config_bak, "w"), allow_unicode = True)             

    def assign_color(self) -> None:
        
        def _get_colors(nodes, name1, name2):
            colors = [(self.config[node].get(name1, None), self.config[node].get(name2, None)) for node in nodes]
            colors = [v for v in colors if all(v)]
            if not colors and self.old_config is not None:
                colors = [(self.old_config[node].get(name1, None), self.old_config[node].get(name2, None)) for node in nodes if node in self.old_config]
                colors = [v for v in colors if all(v)]
            return colors
        
        input_nodes = get_nodes_by_analysis(self.config, "input") | get_nodes_by_analysis(self.config, "tcga_input")
        colors = _get_colors(input_nodes, "control_color", "treat_color")
        if colors:
            control_color, treat_color = colors[0]
        else:
            control_color, treat_color = self.color_assigner.get_group_color()
        for node in input_nodes:
            self.config[node]["control_color"] = control_color
            self.config[node]["treat_color"] = treat_color
                    
            
        low_high_nodes = get_nodes_by_analysis(self.config, "VERIFY") | get_nodes_by_analysis(self.config, "External") | get_nodes_by_analysis(self.config, "PrognosticModel")
        colors = _get_colors(low_high_nodes, "low_color", "high_color")
        if colors:
            low_color, high_color = colors[0]
        else:
            low_color, high_color = self.color_assigner.get_group_color()
        for node in low_high_nodes:
            self.config[node]["high_color"] = high_color
            self.config[node]["low_color"] = low_color
        
        low_high_nodes = get_nodes_by_analysis(self.config, "SSGSEA")
        colors = _get_colors(low_high_nodes, "low_color", "high_color")
        if colors:
            low_color, high_color = colors[0]
        else:
            low_color, high_color = self.color_assigner.get_group_color()
        for node in low_high_nodes:
            self.config[node]["high_color"] = high_color
            self.config[node]["low_color"] = low_color
            
        sub_nodes = get_nodes_by_analysis(self.config, "CONSENSUS")
        colors = _get_colors(sub_nodes, "color1", "color2")
        if colors:
            color1, color2 = colors[0]
        else:
            color1, color2 = self.color_assigner.get_group_color()
        for node in sub_nodes:
            self.config[node]["color1"] = color1
            self.config[node]["color2"] = color2

        aidd_nodes = get_nodes_by_analysis(self.config, "AIDD")
        colors = _get_colors(aidd_nodes, "low_color", "high_color")
        if colors:
            low_color, high_color = colors[0]
        else:
            low_color, high_color = self.color_assigner.get_low_high_color()
        for node in aidd_nodes:
            self.config[node]["low_color"] = low_color
            self.config[node]["high_color"] = high_color

        colors = _get_colors(aidd_nodes, "aidd_mimic_color", "aidd_reverse_color")
        if colors:
            mimic_color, reverse_color = colors[0]
        else:
            mimic_color = self.color_assigner.get_group_color()
            reverse_color = self.color_assigner.get_group_color()
        for node in aidd_nodes:
            self.config[node]["aidd_mimic_color"] = mimic_color
            self.config[node]["aidd_reverse_color"] = reverse_color

        sc_input_nodes = get_nodes_by_analysis(self.config, "sc_input")
        if sc_input_nodes:
            sc_colors = _get_colors(sc_input_nodes, "control_color", "treat_color")
            if sc_colors:
                sc_control_color, sc_treat_color = sc_colors[0]
            elif input_nodes:
                sc_control_color, sc_treat_color = control_color, treat_color
            else:
                sc_control_color, sc_treat_color = self.color_assigner.get_group_color()
            for node in sc_input_nodes:
                self.config[node]["control_color"] = sc_control_color
                self.config[node]["treat_color"] = sc_treat_color

        color_heat = self.config["global"].get("color_heat", None)
        if not color_heat:
            if self.old_config is not None:
                color_heat = self.old_config["global"].get("color_heat", None)
            if not color_heat:
                low_color, high_color = self.color_assigner.get_low_high_color()
                color_heat = [low_color, "white",high_color]
        self.config["global"]["color_heat"] = color_heat
        
        color_panel = self.config["global"].get("color_panel", None)
        if not color_panel:
            if self.old_config is not None:
                color_panel = self.old_config["global"].get("color_panel", None)
            if not color_panel:
                color_panel = self.color_assigner.get_colors()
        self.config["global"]["color_panel"] = color_panel
        
class TcgaCache(Cache):
    def load(self) -> None:
        if not Path(self.cahe_path).exists():
            Path(self.cahe_path).mkdir(parents=True)
        self.dag = parse_node(config=self.config)
        self.old_dag = load_dag(self.dag_path)
        if Path(self.config_bak).exists():
            self.old_config = yaml.load(open(self.config_bak), Loader=yaml.FullLoader)
            input_node = self.old_config.get("tcga_input", None)
            input_change = False
            # if input_node["disease_name"] and input_node["gse"]:
            if input_node:
                if input_node["tcga"]:
                    if self.config["tcga_input"]["tcga"] != input_node["tcga"]:
                        input_change = True
                else:
                    if self.config["tcga_input"]["group"] != input_node["group"]:
                        input_change = True
            if input_change:
                print(f"input 模块的输入数据发生改变，是否删除缓存参数（配色等）和结果（{self.result_path}）？")
            common_change = False
            if self.old_dag is not None:
                common_change =  not compare_dag_graph(self.old_dag, self.dag, print_diff=True)
                if common_change:
                    print(f"缓存流程与配置不一致，是否删除缓存参数（配色等）和结果（{self.result_path}）？")
            if input_change or common_change:
                s = input("输入 y 删除缓存参数，其他任意输入默认保留：")
                if s == "y":
                    self.old_config = None
                    self.old_dag = None
                    if Path(self.result_path).exists():
                        shutil.rmtree(self.result_path)
                    # raise ValueError("缓存流程与配置不兼容，请检查配置文件或者删除缓存文件：{}".format(self.dag_path))
        else:
            self.old_config = None 

    def assign_color(self) -> None:
        
        def _get_colors(nodes, name1, name2):
            colors = [(self.config[node].get(name1, None), self.config[node].get(name2, None)) for node in nodes]
            colors = [v for v in colors if all(v)]
            if not colors and self.old_config is not None:
                colors = [(self.old_config[node].get(name1, None), self.old_config[node].get(name2, None)) for node in nodes if node in self.old_config]
                colors = [v for v in colors if all(v)]
            return colors
        
        input_nodes = get_nodes_by_analysis(self.config, "tcga_input")
        colors = _get_colors(input_nodes, "control_color", "treat_color")
        if colors:
            control_color, treat_color = colors[0]
        else:
            control_color, treat_color = self.color_assigner.get_group_color()
        for node in input_nodes:
            self.config[node]["control_color"] = control_color
            self.config[node]["treat_color"] = treat_color
                    
            
        low_high_nodes = get_nodes_by_analysis(self.config, "VERIFY") | get_nodes_by_analysis(self.config, "External")
        colors = _get_colors(low_high_nodes, "low_color", "high_color")
        if colors:
            low_color, high_color = colors[0]
        else:
            low_color, high_color = self.color_assigner.get_group_color()
        for node in low_high_nodes:
            self.config[node]["high_color"] = high_color
            self.config[node]["low_color"] = low_color
        
        low_high_nodes = get_nodes_by_analysis(self.config, "SSGSEA")
        colors = _get_colors(low_high_nodes, "low_color", "high_color")
        if colors:
            low_color, high_color = colors[0]
        else:
            low_color, high_color = self.color_assigner.get_group_color()
        for node in low_high_nodes:
            self.config[node]["high_color"] = high_color
            self.config[node]["low_color"] = low_color
            
        sub_nodes = get_nodes_by_analysis(self.config, "CONSENSUS")
        colors = _get_colors(sub_nodes, "color1", "color2")
        if colors:
            color1, color2 = colors[0]
        else:
            color1, color2 = self.color_assigner.get_group_color()
        for node in sub_nodes:
            self.config[node]["color1"] = color1
            self.config[node]["color2"] = color2
        color_heat = self.config["global"].get("color_heat", None)
        if not color_heat:
            if self.old_config is not None:
                color_heat = self.old_config["global"].get("color_heat", None)
            if not color_heat:
                low_color, high_color = self.color_assigner.get_low_high_color()
                color_heat = [low_color, "white",high_color]
        self.config["global"]["color_heat"] = color_heat
        
        color_panel = self.config["global"].get("color_panel", None)
        if not color_panel:
            if self.old_config is not None:
                color_panel = self.old_config["global"].get("color_panel", None)
            if not color_panel:
                color_panel = self.color_assigner.get_colors()
        self.config["global"]["color_panel"] = color_panel
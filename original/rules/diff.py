import os
from pathlib import Path
from .rule import MyRule, RuleSet, Parameter


class Diff_check_proportion(MyRule):
    def __init__(self, in_map, incorrect_proportion, out, **kwargs):
        self.script = "diff_check_proportion.py"
        self.args = ["{input.in_map}", "{params.incorrect_proportion}", "{output.out}"]


class Diff(MyRule):
    def __init__(self, in_mat, in_map, out_mat, confirm_file, diff, norm, model, check_done, **kwargs):
        self.script = "diff.R"
        self.args = ["{input.in_mat}", "{input.in_map}", "{output.out_mat}", "{params.confirm_file}", "{params.diff}", "{params.norm}", "{params.model}"]
        self._input["check_done"] = check_done


class Diff_filter(MyRule):
    def __init__(self, confirm_file, in_mat, in_gene, out_mat, rdegs, p_name, p_value, cutoff, logfc_cutoff):
        self.script = "diff_filter.R"
        self.args = ["{input.in_mat}", "{input.in_gene}" if in_gene is not None else "None", "{output.out_mat}",
                     "{output.rdegs}", "{params.confirm_file}",  p_name, p_value, cutoff, logfc_cutoff]

class Diff_venn(MyRule):
    def __init__(self, in_mat, rgs, venn, pheno_abbr, color_panel):
        super().__init__()
        self._input = {"mat": in_mat, "rgs": rgs}
        self._output = {"venn": venn}
        self._params = {}
        self.script = "diff_venn.R"
        self.args = ["{input.mat}", "{input.rgs}", "{output.venn}", pheno_abbr, color_panel]


class Diff_volcano(MyRule):
    def __init__(self, in_mat, in_map, out_mat, volcano, p_name, p_value, logfc_cutoff, volcano_top, volcano_gene):
        super().__init__()
        self._input = {"mat": in_mat, "map": in_map}
        self._output = {"mat": out_mat, "volcano": volcano}
        self._params = {}
        self.script = "diff_volcano.R"
        self.args = ["{input.mat}", "{output.mat}", "{output.volcano}",
                     p_name, p_value, logfc_cutoff, volcano_top, volcano_gene]


class Diff_heatmap(MyRule):
    def __init__(self, in_mat, in_map, rdegs, out_mat, heatmap, confirm_file,  top, color_heat):
        self.script = "diff_heatmap.R"
        self.args = ["{input.in_mat}", "{input.in_map}", "{input.rdegs}",
                     "{output.out_mat}", "{output.heatmap}",  top, color_heat, "{params.confirm_file}"]


class Diff_locate(MyRule):
    def __init__(self, confirm_file, in_mat, locate, plot, tax_id):
        self.script = "diff_locate.R"
        self.args = ["{input.in_mat}", "{input.locate}", "{output.plot}", tax_id, "{params.confirm_file}"]


# class Cp_input(MyRule):
#     def __init__(self, **kwargs):
#         self.script = Path("tools", "cli.py")
#         self.args = ["cp2dir", "{input.in_mat}", "{input.in_map}", "{params.rgs}", "{params.locate}"]
    
#     def post_init(self):
        

class DiffSet(RuleSet):
    def __init__(self, config, name, in_mat: Parameter, in_map: Parameter, in_gene: Parameter, locate,  diff):
        super().__init__(config)
        config = self.config
        self.report_key = ""
        self.ai_file = {
            "a":  Path("assets", "ai", "Z-DiffAnalysis.ai"),
            "b":  Path("assets", "ai", "Z-DiffAnalysis_Chromosome.ai")
        }
        if in_mat.mark != "counts" and diff in ["deseq2", "edgeR"]:
            raise ValueError("{} not support {} matrix.".format(diff, in_mat.mark))
        if in_mat.mark == "counts" and diff == "limma":
            raise ValueError("{} only support counts matrix.".format(diff))
        incorrect_proportion = str(config[name].get("incorrect_proportion", "FALSE"))
        check_done = Path("temp", "proportion_check.done")
        check = Diff_check_proportion(
            in_map=in_map,
            incorrect_proportion=incorrect_proportion,
            out=check_done
        )
        self.append(check)
        diff = Diff(
            in_mat=in_mat,
            in_map=in_map,
            out_mat=Path("temp", "Diffanalysis.csv"),
            confirm_file=self.confirm_file,
            diff=diff,
            norm=config[name]["norm"],
            model=config[name]["model"],
            check_done=check._export["out"]
        )
        self.append(diff)
        self.out_mat= diff._export["out_mat"]

        logfc_cutoff = config[name]["logfc_cutoff"]
        p_name = config[name]["p_name"]
        p_value = config[name]["p_value"]
        diff_filter = Diff_filter(
            in_mat=diff._export["out_mat"], in_gene=in_gene,
            out_mat=Path("temp", "Diffanalysis_logFC={}_{}={}.csv".format(
                logfc_cutoff, p_name, p_value
            )),
            rdegs=Path("temp", "rdegs_logFC={}_{}={}.csv".format(
                logfc_cutoff, p_name, p_value
            )),
            confirm_file=self.confirm_file,
            cutoff=config[name]["cutoff"],
            logfc_cutoff=logfc_cutoff,
            p_name=p_name,
            p_value=p_value
        )
        self.append(diff_filter)
        self.degs = diff_filter._export["out_mat"]
        self.rdegs = diff_filter._export["rdegs"]

        if in_gene is not None:
            diff_venn = Diff_venn(
                diff_filter._export["out_mat"], in_gene,
                venn=Path("output", "2-VennPlot.pdf"),
                pheno_abbr=config[name]["pheno_name"],
                color_panel=",".join(config["global"]["color_panel"])
            )
            self.append(diff_venn)

        diff_volcano = Diff_volcano(
            diff._export["out_mat"], in_map,
            out_mat=Path("output", "1-Vocano.csv"),
            volcano=Path("output", "1-VolcanoPlot.pdf"),
            p_name=p_name,
            p_value=p_value,
            logfc_cutoff=logfc_cutoff,
            volcano_top=config[name]["volcano_top"],
            volcano_gene=config[name]["volcano_gene"]
        )
        self.append(diff_volcano)

        diff_heatmap = Diff_heatmap(
            in_mat=in_mat, in_map=in_map, rdegs=diff_filter._export["rdegs"],
            out_mat=Path("output", "3-Heatmap.csv"),
            heatmap=Path("output", "3-Heatmap.pdf"),
            top=config[name]["top"],
            color_heat=",".join(config["global"]["color_heat"]),
            confirm_file=self.confirm_file
        )
        self.append(diff_heatmap)

        if locate:
            diff_locate = Diff_locate(
                in_mat=diff_filter._export["rdegs"],
                locate=Path(locate),
                plot=Path("output", "4-Chromosome.pdf"),
                tax_id=config["global"]["tax_id"],
                confirm_file=self.confirm_file
            )
            self.append(diff_locate)


class DiffSets(RuleSet):
    def __init__(self, config, name, in_mat, in_map, in_gene, locate):
        super().__init__(config)
        config = self.config
        self.key = "Diffanalysis"
        self.report_key = "Blank"

        diff_method = config[name]["diff"]
        diff = DiffSet(config, name, in_mat, in_map, in_gene, locate,
                       diff=diff_method)
        diff.out_prefix = Path("diff_" + diff_method)
        diff.name = "diff_" + diff_method
        diff.ext_sub_name(suffix="_"+diff_method)
        diff.report_key = "Diffanalysis"
        self.diff = diff
        self.append(diff)
        self.out_mat = diff.out_mat
        self.degs = diff.degs
        self.rdegs = diff.rdegs
        self.out_gene = diff.rdegs

        stat_method = config[name]["stat"]
        self.stat = None
        if stat_method:
            stat = DiffSet(config, name, in_mat, in_map, in_gene, locate,
                           diff=stat_method)
            stat.out_prefix = Path("diff_" + stat_method + "_test")
            stat.name = "diff_" + stat_method
            stat.ext_sub_name(suffix=stat_method)
            self.append(stat)
            self.stat = stat

library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

edf =  io$read_table('signalling_data_log2fc_to_bsa_ctrl.txt', header=TRUE) %>%
	tidyr::gather(measure, value, AKT:Stat3) %>%
    filter(#time == "0.5h",
           treatment != "BSA") %>%
    filter((treatment == "4OHT" & measure == "MEK") |
           (treatment == "IFN" & measure == "Stat3") |
           (treatment == "LY" & measure == "AKT") |
           (treatment == "TGF" & measure == "Smad2") |
           (treatment == "TNF" & measure == "IkBa")) %>%
    group_by(treatment) %>%
    mutate(p.value = t.test(value)$p.value / 2) %>% # one-sided
    mutate(p.value = sprintf("*\np %.2g", p.value)) %>%
    ungroup()

pdf("phospho.pdf", width=6, height=3)

ggplot(edf, aes(x=treatment, y=value, fill=measure)) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    geom_boxplot() +
    geom_point(aes(shape=time), size=3) +
    geom_text(aes(label=p.value), y=-1.5) +
    xlab("Perturbation") +
    ylab("log FC over BSA")

dev.off()

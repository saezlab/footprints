library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

edf =  io$read_table('signalling_data_log2fc_to_bsa_ctrl.txt', header=TRUE) %>%
	tidyr::gather(measure, value, AKT:Stat3) %>%
    filter(treatment != "BSA") %>% #,
#           measure != "cJun") %>%
    group_by(treatment) %>%
    mutate(value = value / max(abs(value))) %>%
    ungroup()

edf2 = edf %>%
    group_by(treatment, measure, time) %>%
    summarize(p.value = t.test(value)$p.value / 2, # one-sided
              value = mean(value),
              label = sprintf("*\np %.2g", p.value)) %>%
    ungroup()

edf3 = edf %>%
    group_by(treatment, measure) %>%
    summarize(p.value = t.test(value)$p.value / 2, # one-sided
              value = mean(value),
              label = sprintf("*\np %.2g", p.value)) %>%
    ungroup()

pdf("hek_phospho.pdf")

edf2 %>%
    mutate(label = ifelse(p.value < 0.05, "*", "")) %>%
    plt$matrix(value ~ treatment + measure, label=label,
               symmetric=TRUE, reverse_color=TRUE) +
    coord_fixed() +
    facet_wrap(~time, ncol=1) +
    ggtitle("* p < 0.05")

edf3 %>%
    mutate(label = ifelse(p.value < 0.05 & abs(value) > 0.3, "*", "")) %>%
    plt$matrix(value ~ treatment + measure, label=label,
               symmetric=TRUE, reverse_color=TRUE) +
    coord_fixed() +
    ggtitle("* p < 0.05 and min 30% max response")

ggplot(edf, aes(x=treatment, y=value, fill=measure)) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    geom_boxplot() +
#    geom_point(aes(shape=time), size=3) +
#    geom_text(aes(label=p.value), y=-1.5) +
    xlab("Perturbation") +
    ylab("log FC over BSA") +
    facet_wrap(~time, ncol=1)

dev.off()

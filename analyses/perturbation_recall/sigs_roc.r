io = import('io')
st = import('stats')
roc = import('./roc')

scoredf = io$load('sigs_scores.RData') %>%
    roc$scores2df() %>%
    mutate(inferred = sub("\\..*$", "", signature))

# per signature, how well do we infer perturbed pathway?
roc = scoredf %>%
    na.omit() %>%
    mutate(matched = perturbed == inferred) %>%
    group_by(inferred, signature) %>%
    do(st$roc(., "score", "matched")) %>%
    ungroup()

width=1
random_line = data.frame(x=c(0,1), y=c(0,1), method=roc$signature[1])
ggplot(roc, aes(x=FPR, y=TPR)) +
	geom_line(aes(x=x, y=y), data=random_line, color="grey", linetype="dashed", size=width) +
	geom_step(size=width, stat="summary", fun.y=median) +
	geom_step(size=width, stat="summary", fun.y=max, linetype="dotted") +
	geom_step(size=width, stat="summary", fun.y=min, linetype="dotted") +
	coord_fixed() +
	facet_wrap(~inferred) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

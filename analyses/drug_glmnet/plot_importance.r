b = import('base')
io = import('io')
ar = import('array')
plt = import('plot')

result = io$load('model_mse.RData')
rmat = ar$construct(mse.test.mean ~ dset + drugs + subset, result, fun.aggregate=function(x) x[1])

rmat = apply(rmat, c(2,3), FUN=function(x) dimnames(rmat)[[1]][which(x==min(x))] %or% NA)

rsum = apply(rmat, 2, FUN=table)

rsum
rowSums(rsum)

pdf("importance.pdf", paper="a4r")
on.exit(dev.off)
df = melt(rmat)
ggplot(df, aes(value, fill=value)) + geom_bar() + ggtitle("all tissues combined")
ggplot(df, aes(value, fill=value)) + geom_bar() + facet_wrap(~Var2) + ggtitle("per tissue")

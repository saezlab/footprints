io = import('io')
lincs = import('data/lincs')

# read index file
index = io$load('index.RData')

# read all relevant expression from LINCS CMap
expr = lincs$get_z(index$distil_id, rid=lincs$bing, map.genes="hgnc_symbol")

# save in object
save(expr, file="expr.RData")

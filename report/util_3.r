io = import('io')
vp = import('../analyses/drug_assocs/plot')

volcano = function(fid) {
    fp = io$file_path('assocs_mapped', fid, ext=".RData")
    assocs = vp$load_fun(fp)$assocs.pan
    vp$plot_pancan(assocs)
}

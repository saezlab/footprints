.dp = import_package_('dplyr')
load('result_call_before_mapply.RData')
load('index_debug.RData')

concat = function(r,n) {
	if (is.null(r))
		list(..id=n)
	else {
		r = as.list(r)
		c(r, ..id=list(rep(n, length(r[[1]]))))
	}
}

result = mapply(concat, result, as.character(seq_along(result)),
				USE.NAMES=TRUE, SIMPLIFY=FALSE) %>%
	.dp$bind_rows() %>%
	.dp$left_join(index, ., by="..id") %>%
	.dp$select(-..id)
result$..id = NULL
colnames(result)[colnames(result) == ""] = "result"


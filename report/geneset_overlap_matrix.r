<<setup_genesets, include=FALSE>>=
overlap = import('../analyses/signature_overview/geneset_overlap')
@

\begin{figure}[H]
<<genesets, echo=FALSE, fig.width=19, fig.height=15>>=
overlap$geneset_overlap_matrix()
@
\caption{Difference in gene sets between Perturbation-response genes and Gene
    Ontology/Reactome. Perturbation-response genes are different to pathway
        members and gene annotations, with between 0 and 4 genes in common.
        However, pathway annotation is quite different from gene annotations
        here as well.}
\label{fig:genesets}
\end{figure}

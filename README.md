## note

snakemake find_hvg --unlock 

script：
    - 不用使用format-string 来生成
    - 不能使用lamada表达式，函数
    - 他的字符串直接可以用大括号

## 调用图

snakemake --dag | dot -Tpdf > dag.pdf

--rulegraph：显示规则之间的依赖关系，而不是具体的文件。
--filegraph：显示文件之间的依赖关系。

snakemake --rulegraph | dot -Tpdf > rulegraph.pdf

snakemake --filegraph | dot -Tpdf > filegraph.pdf ##




snakemake -n
snakemake --export-cwl workflow.cwl
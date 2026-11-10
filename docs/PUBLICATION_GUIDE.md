# Guia para publicação do repositório

## Antes de publicar no GitHub ou arquivar DOI

- Substituir os campos de autoria e referência em `CITATION.cff`.
- Escolher uma licença com todos os autores/detentores de direitos (ex.: MIT) e substituir `LICENSE-NOTICE.txt` por `LICENSE`.
- Revisar `legacy/Metotologia_tese_original.R` e remover caminhos locais ou informações privadas **somente da cópia de publicação**, mantendo uma cópia imutável privada do original.
- Fixar versões de R e dependências (por exemplo, com `renv`) e executar testes em sistemas suportados.
- Disponibilizar, se permitido, pelo menos uma cena de demonstração georreferenciada, seu manifesto de partições, a exportação limpa do FracLac e um resultado de referência.
- Validar numericamente os descritores e os rasters com saídas originais utilizadas no artigo.
- Tornar explícita a diferença entre a configuração **original** `criterion='SSErr'` e o default **novo** `criterion='RMSE_CV'`.
- Documentar os q amostrados, valores limites, interpolação, parâmetros do variograma e critérios de recorte.

## Redação sugerida para a seção de disponibilidade de código (adaptar antes do depósito)

> The R scripts implementing partition indexing, geospatial attribution of FracLac-derived multifractal descriptors, empirical semivariogram fitting, and ordinary kriging are available in the accompanying repository [URL/DOI to be supplied]. Multifractal spectra were computed externally with FracLac in ImageJ. The repository includes reproducibility checks and an example workflow; original satellite scenes and selected FracLac outputs are [availability conditions to be supplied].

## Posicionamento científico

O repositório oferece um **workflow reproduzível de pós-processamento e espacialização** de métricas obtidas com FracLac, e não um novo algoritmo de FracLac, uma demonstração de equivalência com a implementação Java ou um método de krigagem inédito.

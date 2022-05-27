# RGR-bio-review

Code to retrive biological records from GBIF and OBIS databases in the Rio Grande Rise, used in the paper:

Montserrat, F., Guilhon, M., Corrêa, P. V. F., et al. (2019). Deep-sea mining on the Rio Grande Rise (Southwestern Atlantic): A review on environmental baseline, ecosystem services and potential impacts. Deep Sea Research Part I: Oceanographic Research Papers, 145, 31–58. https://doi.org/10.1016/j.dsr.2018.12.007.

`table_get.R` gets the original data retrived from the databases *obis_bak.csv.gz* and *0101012-160910150852091.zip*, clean it, and save as *obis.csv.gz* and *gbif.csv.gz*<br>
`table_generate.R` joins data from both databases and filters only animal records

If you use this code, please cite our paper.

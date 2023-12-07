# CHEST-CRITCARE-D-23-00091

Publicly-available code for Dr. Ellen Burnham's manuscript, titled: "Prevalence of Alcohol Use Characterized by Phosphotidylethanol and Its Impact on Outcomes in Hospitalized Patients with Respiratory Failure Before and During the COVID-19 Pandemic."

Coauthors: Ellen Burnham, Grace Perry, Patrick Offner, Ryan Ormesher, Jeanette Gaydos, Suzanne Slaughter, Carrie Higgins, Jeffrey McKeean, Raymond Pomponio, Ryan Peterson, Sarah Jolly

## Code

All code used to analyze the data and produce the results in the paper can be found under `Scripts/`. `R`` version 4.2.1 (2022-06-23) was used.

## Figures

The figures produced by `R` can be found under `Figures/`.

## Output

The file `Output/regenerate_all.qmd` is used to aggregate all of the tables and figures used in the paper. It is a Quarto [markdown](https://quarto.org/docs/computations/r.html) file that runs all of the `R` code sequentially and gathers output in an *html* document.

The resulting file, `Output/regenerate_all.html`, can be opened in any web browser.
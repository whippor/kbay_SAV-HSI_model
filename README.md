# Kachemak Bay Submerged Aquatic Vegetation Habitat Suitability Index Model

###### README LAST UPDATED: 2025/04/29

------------------------------------------------------------------------

## Overview

------------------------------------------------------------------------

This is a repository under active development to model habitat suitability indices (HSI) for submerged aquatic vegetation (SAV) including canopy kelp, understorey kelp, and seagrass in Kachemak Bay, Alaska.

### Objectives

-   Create a validated habitat suitability index model for major SAV categories in Kachemak Bay.

-   Provide publicly available data to identify critical SAV types and potential locations.

### Guiding Questions

-   What is the theoretical total area that canopy kelps, understorey kelps, and seagrass could colonize in the bay?

-   How do HSI values align with past and present observations of seagrass and kelp?

### Current Status

This project is in active development and will change as additional datasets and analyses are included. The primary short-term goals for the project include:

-   ~~Streamline the code volume by reducing duplicated files created for each SAV type to a single reference file (initial code architecture was highly redundant to allow modularity of component SAV parts).~~ **Completed: 2025/03/06**

-   Identify additional environmental parameters and datasets relevant to habitat suitability, especially oceanographic parameters including temperature, salinity, and nutrients.

-   Ensure proper functioning of local fetch calculation code.

-   ~~Validate all TAM table values with references.~~ **Completed: 2025/04/28**

-   Perform model selection and sensitivity analyses on models and submodels to identify the number and identity of parameters that explain the most model variance.

-   Conduct model validation through field survey methods across all SAV types, HSI values, and locales.

-   Create markdown summary document with major methods, inputs, and outputs.

------------------------------------------------------------------------

## Running The Code

------------------------------------------------------------------------

Once cloned, the code can by run piecemeal for each SAV group and parameter, or the combined HSI model for each SAV type can be run with the **\*\_HSI_model.R** code. Or you can run all data prep and models using the **\*\_HSI_runall.R** code blocks for each SAV. The HSI model for all SAV combined can be run with **721_combined_HSI_model.R.** Alternatively, you can run all the models and submodels using the **731_HSI_runall.R** code block.

| File Name Prefix Value | Content                                     |
|------------------------|---------------------------------------------|
| 000                    | functions required for the analysis         |
| 00\*                   | base layer import, prep, and/or calculation |
| 1\*\*                  | understorey kelp code components            |
| 2\*\*                  | canopy kelp code components                 |
| 3\*\*                  | seagrass code components                    |
| 7\*\*                  | multi-model execution                       |
| 888                    | additional analysis and summary             |
| 999                    | scratch pad                                 |
| U\*\*                  | uncertainty analysis components             |
| V\*\*                  | validation analysis components              |
| \*1\*                  | submodel                                    |
| \*2\*                  | SAV-specific HSI model                      |
| \*31                   | "run all" execution                         |
| \*11                   | bathymetry submodel                         |
| \*12                   | substrate submodel                          |
| \*13                   | fetch submodel                              |
| \*14                   | velocity submodel                           |
| \*15                   | canopy presence submodel                    |
| \*16                   | seagrass presence submodel                  |
| \*21                   | SAV-specific HSI model                      |

: Code type by naming convention.

Currently, the initial fetch calculation code is not operational due to an incompatibility between the fetchR package and the current version of R. However, the original fetch map has been placed in each relevant data directory. In this way the model can still be run for each type of SAV, but fetch will not be calculated locally.

The gitignore file can be changed to ignore derived model data including .tifs if desired.

------------------------------------------------------------------------

## Methods

------------------------------------------------------------------------

### Data Collection

Existing datasets for Kachemak Bay were identified for the current model and synthesized to provide bathymetry, substrate type, and SAV presence (canopy kelps, seagrass) layers.

### Bathymetry

Bathymetry data were generated by Martin Renner from: NOAA’s National Geophysical Data Center Kachemak Bay DEM for tsunami modeling (NOAA 2010); NOAA essential fish habitat smooth sheets (curated by Mark Zimmerman, NMFS), and GMRTv4.2 (2024). Grids were produced by reprojecting and bi-cubically up-sampling each existing grid where appropriate, and manually filling in all missing values. Some irregularities and artifacts are expected at the boundaries between coverages. Visual inspections were performed to address some issues, but others likely remain.

### Substrate Type

Substrate type data were taken from NCEI archives as a thematic map of benthic habitat classification in Kachemak Bay created by multibeam echosounder data, automated image processing techniques, and on-screen editing (NOAA 2020).

### Fetch

Fetch was calculated from the full bathymetry layer as reference for land and water features using functions modified from the 'get_fetch()' command in the fetchR package. See the [fetchR](http://cran.nexr.com/web/packages/fetchR/index.html) documentation for details.

### SAV Presence

Both kelp and seagrass presence data were extracted from the AOOS Workspace and included kelp cover polygons of Kachemak Bay for 2000-2002, and seagrass polygons from 2005.

------------------------------------------------------------------------

## References

------------------------------------------------------------------------

NOAA National Geophysical Data Center. 2010: Kachemak Bay, Alaska 1/3 arc-second MHW Coastal Digital Elevation Model. NOAA National Centers for Environmental Information. Accessed [date].

NCCOS Mapping: Seafloor Mapping Products for Kachemak Bay, Cook Inlet, AK, from 2005-07-06 to 2017-07-19. 2020. NOAA National Centers for Environmental Information.

R Core Team (2024). \_R: A Language and Environment for Statistical Computing\_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

Hijmans R (2024). \_terra: Spatial Data Analysis\_. R package version 1.7-78, <https://CRAN.R-project.org/package=terra>.

Madden CJ, Grossman DH, Goodin KL, Dethier MN. 2005. Coastal and marine systems of North America: framework for an ecological classification standard: version II.

Carrillo, C., S.K. McKay, and T.S. Swannack. 2020. Ecological Model Development: Toolkit for interActive Modeling (TAM). ERDC/TN EMRRP-SR-90. Vicksburg, MS: US Army Engineer Research and Development Center.

------------------------------------------------------------------------

### TAM References

| SAV              | Parameter      | Citation                                                                                                                                                                                                                                                                   |
|------------------|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canopy kelp      | depth          | Springer YP, Hays CG, Carr MH, Mackey MR. Toward ecosystem-based management of marine macroalgae—The bull kelp, *Nereocystis luetkeana*. Oceanography and marine biology. 2010 May 12;48:1.                                                                                |
| canopy kelp      | substrate      | Dayton PK. Ecology of kelp communities. Annual review of ecology and systematics. 1985 Jan 1:215-45.                                                                                                                                                                       |
| canopy kelp      | fetch          | Smale DA, Burrows MT, Evans AJ, King N, Sayer MD, Yunnie AL, Moore PJ. Linking environmental variables with regional-scale variability in ecological structure and standing stock of carbon within UK kelp forests. Marine Ecology Progress Series. 2016 Jan 19;542:79-95. |
| canopy kelp      | temperature \* | Weigel BL, Small SL, Berry HD, Dethier MN. Effects of temperature and nutrients on microscopic stages of the bull kelp (*Nereocystis luetkeana*, Phaeophyceae). Journal of phycology. 2023 Oct;59(5):893-907.                                                              |
| canopy kelp      | salinity \*    |                                                                                                                                                                                                                                                                            |
| understorey kelp | depth          | Bekkby T, Smit C, Gundersen H, Rinde E, Steen H, Tveiten L, Gitmark JK, Fredriksen S, Albretsen J, Christie H. The abundance of kelp is modified by the combined impact of depth, waves and currents. Frontiers in Marine Science. 2019 Aug 6;6:475.                       |
| understorey kelp | substrate      | Dayton PK. Ecology of kelp communities. Annual review of ecology and systematics. 1985 Jan 1:215-45.                                                                                                                                                                       |
| understorey kelp | fetch          | Bekkby T, Moy FE. Developing spatial models of sugar kelp (*Saccharina latissima*) potential distribution under natural conditions and areas of its disappearance in Skagerrak. Estuarine, Coastal and Shelf Science. 2011 Dec 20;95(4):477-83.                            |
| understorey kelp | temperature \* | Bolton JJ, Lüning K. Optimal growth and maximal survival temperatures of Atlantic Laminaria species (Phaeophyta) in culture. Marine Biology. 1982 Jan;66:89-94.                                                                                                            |
| understorey kelp | salinity \*    | Spurkland T, Iken K. Salinity and irradiance effects on growth and maximum photosynthetic quantum yield in subarctic *Saccharina latissima* (Laminariales, Laminariaceae).                                                                                                 |
| seagrass         | depth          | Thom RM, Southard SL, Borde AB, Stoltz P. Light requirements for growth and survival of eelgrass (*Zostera marina* L.) in Pacific Northwest (USA) estuaries. Estuaries and Coasts. 2008 Nov;31:969-80.                                                                     |
| seagrass         | substrate      | Larkum AW, Orth RJ, Duarte CM. Seagrasses: biology, ecology and conservation. Phycologia. 2006;45(5):5.                                                                                                                                                                    |
| seagrass         | fetch          | Oreska MP, McGlathery KJ, Wiberg PL, Orth RJ, Wilcox DJ. Defining the *Zostera marina* (eelgrass) niche from long-term success of restored and naturally colonized meadows: Implications for seagrass restoration. Estuaries and coasts. 2021 Mar;44(2):396-411.           |
| seagrass         | temperature \* | Nejrup LB, Pedersen MF. Effects of salinity and water temperature on the ecological performance of *Zostera marina*. Aquatic Botany. 2008 Apr 1;88(3):239-46.                                                                                                              |
| seagrass         | salinity \*    | Nejrup LB, Pedersen MF. Effects of salinity and water temperature on the ecological performance of *Zostera marina*. Aquatic Botany. 2008 Apr 1;88(3):239-46.                                                                                                              |

\* Not currently included in HSI calculations



# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
#install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",16.463,4.4584,0.995,0.000,"immune system process"),
c("GO:0007169","transmembrane receptor protein tyrosine kinase signaling pathway",3.976,10.3054,0.870,0.000,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0044700","single organism signaling",36.544,10.3809,0.926,0.301,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:1901698","response to nitrogen compound",5.326,6.7144,0.926,0.458,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0023051","regulation of signaling",17.946,9.5376,0.891,0.468,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0051716","cellular response to stimulus",40.358,7.8665,0.932,0.294,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0048017","inositol lipid-mediated signaling",1.212,8.0888,0.880,0.349,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0048583","regulation of response to stimulus",21.610,10.5952,0.882,0.183,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0035556","intracellular signal transduction",15.101,8.1798,0.841,0.440,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0010646","regulation of cell communication",17.651,9.8356,0.887,0.256,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0009628","response to abiotic stimulus",6.705,3.8327,0.953,0.141,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0007267","cell-cell signaling",9.025,4.3768,0.923,0.373,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0071363","cellular response to growth factor stimulus",3.618,7.4828,0.905,0.165,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0099531","presynaptic process involved in chemical synaptic transmission",0.894,7.3161,0.927,0.161,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0070848","response to growth factor",3.762,8.2411,0.919,0.456,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0009611","response to wounding",3.664,4.7620,0.954,0.110,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0009605","response to external stimulus",12.043,5.3778,0.949,0.197,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0007166","cell surface receptor signaling pathway",15.868,9.7852,0.840,0.447,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0042221","response to chemical",23.993,7.2255,0.942,0.307,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0010033","response to organic substance",16.515,8.8210,0.912,0.424,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0014066","regulation of phosphatidylinositol 3-kinase signaling",0.837,9.5817,0.866,0.181,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0009719","response to endogenous stimulus",9.175,12.6234,0.951,0.128,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0033554","cellular response to stress",10.496,5.3279,0.931,0.438,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0031667","response to nutrient levels",2.383,4.2069,0.957,0.103,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0006950","response to stress",21.310,9.7328,0.943,0.233,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0010243","response to organonitrogen compound",4.697,7.0443,0.888,0.474,"transmembrane receptor protein tyrosine kinase signaling pathway"),
c("GO:0009987","cellular process",90.329,4.5834,0.999,0.000,"cellular process"),
c("GO:0023052","signaling",36.613,10.3497,0.996,0.000,"signaling"),
c("GO:0032501","multicellular organismal process",41.143,8.6003,0.997,0.000,"multicellular organismal process"),
c("GO:0032502","developmental process",33.982,9.6144,0.996,0.000,"developmental process"),
c("GO:0040007","growth",5.447,3.8697,0.995,0.000,"growth"),
c("GO:0040011","locomotion",9.452,5.1574,0.995,0.000,"locomotion"),
c("GO:0044699","single-organism process",78.194,9.0531,0.999,0.000,"single-organism process"),
c("GO:0050896","response to stimulus",49.302,8.6925,0.997,0.000,"response to stimulus"),
c("GO:0051179","localization",36.018,10.2976,0.996,0.000,"localization"),
c("GO:0065007","biological regulation",67.069,8.5436,0.998,0.000,"biological regulation"),
c("GO:0071840","cellular component organization or biogenesis",36.532,6.6968,0.996,0.000,"cellular component organization or biogenesis"),
c("GO:0099003","vesicle-mediated transport in synapse",0.767,10.9245,0.906,0.000,"vesicle-mediated transport in synapse"),
c("GO:0006810","transport",29.215,8.8962,0.878,0.496,"vesicle-mediated transport in synapse"),
c("GO:0006836","neurotransmitter transport",1.114,6.7282,0.897,0.193,"vesicle-mediated transport in synapse"),
c("GO:0032879","regulation of localization",14.409,11.4461,0.851,0.211,"vesicle-mediated transport in synapse"),
c("GO:0006898","receptor-mediated endocytosis",2.014,7.5986,0.894,0.439,"vesicle-mediated transport in synapse"),
c("GO:0016192","vesicle-mediated transport",11.385,10.0391,0.894,0.318,"vesicle-mediated transport in synapse"),
c("GO:0099504","synaptic vesicle cycle",0.716,10.0964,0.907,0.409,"vesicle-mediated transport in synapse"),
c("GO:0048278","vesicle docking",0.375,5.1096,0.896,0.368,"vesicle-mediated transport in synapse"),
c("GO:0033036","macromolecule localization",16.624,3.9747,0.905,0.385,"vesicle-mediated transport in synapse"),
c("GO:1902578","single-organism localization",20.162,10.0535,0.887,0.369,"vesicle-mediated transport in synapse"),
c("GO:0051674","localization of cell",8.257,4.8386,0.915,0.295,"vesicle-mediated transport in synapse"),
c("GO:0051641","cellular localization",14.905,5.5768,0.906,0.373,"vesicle-mediated transport in synapse"),
c("GO:0051640","organelle localization",2.949,7.2815,0.904,0.417,"vesicle-mediated transport in synapse"),
c("GO:0051650","establishment of vesicle localization",1.431,10.8962,0.879,0.185,"vesicle-mediated transport in synapse"),
c("GO:0048284","organelle fusion",1.321,6.6840,0.972,0.003,"organelle fusion"),
c("GO:0071825","protein-lipid complex subunit organization",0.219,4.0255,0.981,0.127,"organelle fusion"),
c("GO:0016043","cellular component organization",35.684,7.1427,0.968,0.242,"organelle fusion"),
c("GO:0016050","vesicle organization",2.256,5.2924,0.971,0.325,"organelle fusion"),
c("GO:1902589","single-organism organelle organization",10.387,5.1029,0.935,0.436,"organelle fusion"),
c("GO:0061025","membrane fusion",1.425,5.2644,0.964,0.155,"organelle fusion"),
c("GO:0061024","membrane organization",7.299,4.6968,0.973,0.350,"organelle fusion"),
c("GO:0006996","organelle organization",19.411,3.9318,0.969,0.479,"organelle fusion"),
c("GO:0051128","regulation of cellular component organization",13.335,6.0039,0.902,0.420,"organelle fusion"),
c("GO:0006914","autophagy",2.643,4.7496,0.991,0.003,"autophagy"),
c("GO:0007154","cell communication",36.705,11.7545,0.991,0.005,"cell communication"),
c("GO:0022406","membrane docking",1.004,4.5498,0.968,0.014,"membrane docking"),
c("GO:0006066","alcohol metabolic process",1.702,7.3925,0.892,0.015,"alcohol metabolism"),
c("GO:0044711","single-organism biosynthetic process",8.540,4.8153,0.866,0.425,"alcohol metabolism"),
c("GO:0019216","regulation of lipid metabolic process",1.685,6.3307,0.823,0.260,"alcohol metabolism"),
c("GO:0030214","hyaluronan catabolic process",0.087,3.7258,0.945,0.380,"alcohol metabolism"),
c("GO:0006629","lipid metabolic process",7.963,6.5003,0.886,0.321,"alcohol metabolism"),
c("GO:0044281","small molecule metabolic process",12.072,4.7773,0.883,0.464,"alcohol metabolism"),
c("GO:1903578","regulation of ATP metabolic process",0.317,3.8761,0.804,0.378,"alcohol metabolism"),
c("GO:0048514","blood vessel morphogenesis",2.856,7.2495,0.899,0.016,"blood vessel morphogenesis"),
c("GO:0044707","single-multicellular organism process",35.118,9.8508,0.945,0.227,"blood vessel morphogenesis"),
c("GO:0072359","circulatory system development",5.361,7.1051,0.922,0.312,"blood vessel morphogenesis"),
c("GO:0090077","foam cell differentiation",0.196,4.0255,0.943,0.149,"blood vessel morphogenesis"),
c("GO:0051239","regulation of multicellular organismal process",15.268,9.5935,0.904,0.357,"blood vessel morphogenesis"),
c("GO:0044767","single-organism developmental process",33.474,9.6459,0.913,0.454,"blood vessel morphogenesis"),
c("GO:0050793","regulation of developmental process",12.949,11.3635,0.882,0.254,"blood vessel morphogenesis"),
c("GO:0009653","anatomical structure morphogenesis",13.872,6.1537,0.938,0.464,"blood vessel morphogenesis"),
c("GO:0072001","renal system development",1.621,4.5346,0.927,0.263,"blood vessel morphogenesis"),
c("GO:0001655","urogenital system development",1.823,4.3556,0.933,0.267,"blood vessel morphogenesis"),
c("GO:0044763","single-organism cellular process",62.712,12.5376,0.969,0.029,"single-organism cellular process"),
c("GO:0001505","regulation of neurotransmitter levels",1.102,6.5952,0.941,0.034,"regulation of neurotransmitter levels"),
c("GO:0042592","homeostatic process",9.371,5.1864,0.922,0.378,"regulation of neurotransmitter levels"),
c("GO:0042593","glucose homeostasis",1.252,5.1355,0.935,0.287,"regulation of neurotransmitter levels"),
c("GO:2000505","regulation of energy homeostasis",0.110,3.7258,0.916,0.383,"regulation of neurotransmitter levels"),
c("GO:0006091","generation of precursor metabolites and energy",2.129,3.9469,0.959,0.038,"generation of precursor metabolites and energy"),
c("GO:0008283","cell proliferation",11.321,8.5100,0.978,0.041,"cell proliferation"),
c("GO:1901615","organic hydroxy compound metabolic process",2.689,6.3526,0.964,0.045,"organic hydroxy compound metabolism"),
c("GO:0010675","regulation of cellular carbohydrate metabolic process",0.756,7.7375,0.828,0.049,"regulation of cellular carbohydrate metabolism"),
c("GO:0051246","regulation of protein metabolic process",14.651,5.2950,0.825,0.439,"regulation of cellular carbohydrate metabolism"),
c("GO:0043467","regulation of generation of precursor metabolites and energy",0.502,5.7399,0.872,0.141,"regulation of cellular carbohydrate metabolism"),
c("GO:1902893","regulation of pri-miRNA transcription from RNA polymerase II promoter",0.167,4.8041,0.876,0.129,"regulation of cellular carbohydrate metabolism"),
c("GO:0061614","pri-miRNA transcription from RNA polymerase II promoter",0.179,4.6478,0.930,0.292,"regulation of cellular carbohydrate metabolism"),
c("GO:0060255","regulation of macromolecule metabolic process",33.630,5.1343,0.805,0.288,"regulation of cellular carbohydrate metabolism"),
c("GO:0009894","regulation of catabolic process",3.006,4.4353,0.863,0.158,"regulation of cellular carbohydrate metabolism"),
c("GO:0043112","receptor metabolic process",0.917,3.7375,0.952,0.049,"receptor metabolism"),
c("GO:0010506","regulation of autophagy",1.575,3.9547,0.942,0.054,"regulation of autophagy"),
c("GO:0065008","regulation of biological quality",20.202,12.7305,0.937,0.055,"regulation of biological quality"),
c("GO:0008284","positive regulation of cell proliferation",4.859,9.8794,0.838,0.064,"positive regulation of cell proliferation"),
c("GO:0048522","positive regulation of cellular process",27.732,10.7190,0.819,0.376,"positive regulation of cell proliferation"),
c("GO:0030335","positive regulation of cell migration",2.308,6.4724,0.760,0.273,"positive regulation of cell proliferation"),
c("GO:0031328","positive regulation of cellular biosynthetic process",10.167,10.9136,0.736,0.453,"positive regulation of cell proliferation"),
c("GO:0051240","positive regulation of multicellular organismal process",8.327,8.1421,0.821,0.429,"positive regulation of cell proliferation"),
c("GO:0006357","regulation of transcription from RNA polymerase II promoter",10.900,7.7825,0.795,0.463,"positive regulation of cell proliferation"),
c("GO:0045927","positive regulation of growth",1.414,4.6180,0.887,0.222,"positive regulation of cell proliferation"),
c("GO:0045597","positive regulation of cell differentiation",4.893,5.6289,0.793,0.307,"positive regulation of cell proliferation"),
c("GO:0010885","regulation of cholesterol storage",0.075,4.0353,0.880,0.302,"positive regulation of cell proliferation"),
c("GO:0046324","regulation of glucose import",0.369,5.2204,0.842,0.354,"positive regulation of cell proliferation"),
c("GO:0051130","positive regulation of cellular component organization",6.861,6.2299,0.836,0.408,"positive regulation of cell proliferation"),
c("GO:0010743","regulation of macrophage derived foam cell differentiation",0.167,4.3072,0.900,0.474,"positive regulation of cell proliferation"),
c("GO:0043066","negative regulation of apoptotic process",4.749,8.8268,0.876,0.066,"negative regulation of apoptotic process"),
c("GO:0048523","negative regulation of cellular process",24.951,6.6968,0.899,0.434,"negative regulation of apoptotic process"),
c("GO:0048585","negative regulation of response to stimulus",8.003,3.9872,0.873,0.491,"negative regulation of apoptotic process"),
c("GO:0008219","cell death",11.881,9.0931,0.957,0.067,"cell death"),
c("GO:0050790","regulation of catalytic activity",13.993,6.0477,0.920,0.074,"regulation of catalytic activity"),
c("GO:0006793","phosphorus metabolic process",18.702,10.9101,0.943,0.076,"phosphorus metabolism"),
c("GO:0065009","regulation of molecular function",17.288,7.3468,0.939,0.078,"regulation of molecular function"),
c("GO:0009056","catabolic process",11.471,3.9208,0.959,0.079,"catabolism"),
c("GO:0006928","movement of cell or subcellular component",10.987,4.5498,0.958,0.080,"movement of cell or subcellular component"),
c("GO:0048518","positive regulation of biological process",30.969,10.5834,0.916,0.095,"positive regulation of biological process"),
c("GO:0048519","negative regulation of biological process",26.855,6.3063,0.919,0.134,"positive regulation of biological process"),
c("GO:0019222","regulation of metabolic process",35.730,4.7825,0.900,0.177,"positive regulation of biological process"),
c("GO:0050794","regulation of cellular process",60.427,7.0915,0.899,0.274,"positive regulation of biological process"),
c("GO:0050789","regulation of biological process",63.456,8.1752,0.918,0.162,"positive regulation of biological process"),
c("GO:0044710","single-organism metabolic process",24.524,8.1512,0.938,0.099,"single-organism metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()

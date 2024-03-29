% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_species.r
\name{compare_spp_names}
\alias{compare_spp_names}
\title{Compare species names between the sample dataset and the comparison
pusedo-quadrat dataset}
\usage{
compare_spp_names(samp_df, comp_df = NULL, spp_col = "species", verbose = TRUE)
}
\arguments{
\item{samp_df}{A \code{data.frame} containing the sample data that you want to
assign to NVC via comparison with pseudo-quadrat data}

\item{comp_df}{A \code{data.frame} containing the comparison data samples you
want to match your focal samples in \code{samp_df} to. The default value of
\code{NULL} means that the comparison will use the NVC pseudo-quadrat data
included with the package}

\item{spp_col}{A character vector specifying the name of the column in
\code{samp_df} and \code{comp_df}that contain the species names or codes that you
want to match on. The default is for the column to be named \code{species}.}

\item{verbose}{A logical variable determining whether the function should
output text to the display}
}
\value{
A \code{character} vector containing the names of species from the sample
data that were not matched to names in the pseudo-quadrant data
}
\description{
Compare the species names used in the sample dataset against
those used in the pseudo-quadrat dataset and highlight any species names
from the sample not present in the pseudo-quadrat data it is to be compared
against. Mismatched names may be the result of species in the focal data
that are not included in or relevant to the NVC psueod-quadrat, in which
case they can be safely ignored. However it is also possible that they may
be due to differences in nomenclature used in the two datasets and that the
species could be in both but under different names. Ideally the datasets
should both use a  standardised nomenclature, however this may not often be
possible.
}
\details{
Currently the NVC pseudo-quadrat data uses relatively old
nomenclature so there are likely to be issues with mismatched names when
compared against datasets using relatively recent nomenclature. In future
we plant to update the nomenclature used in the nvc psudo-quadrats to a
more modern nomenclature and would like to ideally include or link to
packages that can help with the standardisation of taxnomomies between the
datasets but for now you should try to ensure that where possible your
nomenclature matches the NVC nomenclature as fully as possible.
}
\examples{

#' # Create a test dataset
test_samp = data.frame(ID = 1, species = c("Alchemilla alpina","Anthoxanthum odoratum","Blechnum spicant","Deschampsia flexuosa","Galium saxatile","Luzula sylvatica","Melampyrum pratense","Nardus stricta","Potentilla erecta","Rumex acetosella","Vaccinium myrtillus","Dicranum majus","Hylocomium splendens","Isothecium myosuroides","Plagiothecium undulatum","Pleurozium schreberi","Polytrichum commune","Racomitrium lanuginosum","Rhytidiadelphus loreus","Calypogeia fissa","Unidentified plant","Cladonia bellidiflora","Cladonia coccifera","Cladonia uncialis"))

# Now run names comparison function on this test data
unmatched_spp = compare_spp_names(test_samp)

}

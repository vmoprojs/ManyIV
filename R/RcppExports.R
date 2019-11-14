arMM <- function(A,B)
{
    .Call('_ManyIV_arMM',PACKAGE = 'ManyIV',A,B)
}


eiMM <- function(A,B)
{
    .Call('_ManyIV_eiMM',PACKAGE = 'ManyIV',A,B)
}


eiMu <- function(A,B)
{
    .Call('_ManyIV_eiMMM',PACKAGE = 'ManyIV',A,B)
}

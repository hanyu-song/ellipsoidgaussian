## Resubmission
This is a resubmission. In this version I have:

* I modified the examples so that now only ellipsoid_gaussian potentially runs 
for more than 5 seconds. The function ellipsoid_gasusian needs more time to run 
because the distribution it is sampling from is very complicated - it involves 
numerically approximating the normalising constant of a Fisher-Bingham 
distribution at every iteration for all the data points in the mini-batch, which 
makes the function slow, even though I already coded the slow parts in C++. As a 
result, I put this function in a donotrun block.
* I understand that CRAN checks may identify other issues, including
  - Possibly misspelled words in DESCRIPTION:
        Dunson (8:28)
        Mises (17:5)
        von (16:61)
  - Found the following (possibly) invalid URLs: 
  URL: http://www.jstor.org/stable/20441200
  From: inst/doc/my-vignette.html
  Status: 403
  Message: Forbidden
  
  Below are my explanations:
  - These words are not misspelled. Dunson is the last name of a co-author on
  the cited paper in DESCRIPTION. "von" and "Mises" are part of the name of a
  distribution on the sphere; see
  https://en.wikipedia.org/wiki/Von_Mises–Fisher_distribution for more details.
  - The URL is part of a reference in the vignette. Since I used the method in
  that paper, I thought it should be cited.

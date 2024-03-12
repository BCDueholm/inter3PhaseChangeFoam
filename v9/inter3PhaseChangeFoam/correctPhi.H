correctUphiBCs(U, phi, true);

CorrectPhi
(
    phi,
    U,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU(),
    pressureReference,
    pimple
);

#include "continuityErrs.H"

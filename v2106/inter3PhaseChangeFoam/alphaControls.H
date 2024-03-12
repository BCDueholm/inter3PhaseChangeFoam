const dictionary& alphaControls = mesh.solverDict(alpha1.name());

const label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

const label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

const bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

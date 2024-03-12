const dictionary& alphaControls = mesh.solverDict(alpha1.name());

const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

const bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));


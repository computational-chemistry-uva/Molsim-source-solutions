#include <pybind11/pybind11.h>

#include "hardDisks.h"
#include "mc.h"
#include "md.h"
#include "sample.h"

PYBIND11_MODULE(molsim, m)
{
  pybind11::class_<HardDisks>(m, "HardDisks")
      .def(pybind11::init<size_t, size_t, size_t, double, size_t, double, size_t, bool>(),
           pybind11::arg("numberOfInitCycles"), pybind11::arg("numberOfProdCycles"), pybind11::arg("numberOfParticles"),
           pybind11::arg("maxDisplacement"), pybind11::arg("sampleFrequency"), pybind11::arg("boxSize"),
           pybind11::arg("rdfBins"), pybind11::arg("runStatic"))
      .def("run", &HardDisks::run)
      .def("getRDF", &HardDisks::getRDF)
      .def_readonly("acceptanceRatio", &HardDisks::acceptanceRatio);

  pybind11::class_<MonteCarlo>(m, "MonteCarlo")
      .def(pybind11::init<size_t, size_t, size_t, double, double, double, double, double, double, double, double,
                          size_t, double, double, size_t, size_t>(),
           pybind11::arg("numberOfParticles"), pybind11::arg("numberOfInitCycles"), pybind11::arg("numberOfProdCycles"),
           pybind11::arg("temperature"), pybind11::arg("boxSize"), pybind11::arg("maxDisplacement"),
           pybind11::arg("translationProbability") = 0.0, pybind11::arg("pressure") = 0.0,
           pybind11::arg("volumeProbability") = 0.0, pybind11::arg("maxVolumeChange") = 1.0,
           pybind11::arg("swapProbability") = 0.0, pybind11::arg("sampleFrequency") = 100, pybind11::arg("sigma") = 1.0,
           pybind11::arg("epsilon") = 1.0, pybind11::arg("logLevel") = 0, pybind11::arg("seed") = 12)
      .def("__repr__", &MonteCarlo::repr)
      .def("run", &MonteCarlo::run)
      .def_readonly("pressures", &MonteCarlo::pressures);

  pybind11::class_<MolecularDynamics>(m, "MolecularDynamics")
      .def(pybind11::init<size_t, double, double, double, size_t, size_t, size_t, bool>(),
           pybind11::arg("numberOfParticles"), pybind11::arg("temperature"), pybind11::arg("dt"),
           pybind11::arg("boxSize"), pybind11::arg("sampleFrequency") = 100, pybind11::arg("logLevel") = 0,
           pybind11::arg("seed") = 12, pybind11::arg("useNoseHoover") = false)
      .def_readonly("rdfSampler", &MolecularDynamics::rdfSampler)
      .def_readonly("msdSampler", &MolecularDynamics::msdSampler)
      .def("__repr__", &MolecularDynamics::repr)
      .def("run", &MolecularDynamics::run, pybind11::arg("numberOfSteps"), pybind11::arg("equilibrate"),
           pybind11::arg("outputPDB") = true);
  pybind11::class_<SampleRDF>(m, "SampleRDF").def("getResults", &SampleRDF::getResults);
  pybind11::class_<SampleMSD>(m, "SampleMSD").def("getResults", &SampleMSD::getResults);
}

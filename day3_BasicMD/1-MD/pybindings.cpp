#include <pybind11/pybind11.h>

#include "md.h"
#include "sample.h"

PYBIND11_MODULE(md, m)
{
  pybind11::class_<MolecularDynamics>(m, "MolecularDynamics")
      .def(pybind11::init<size_t, double, double, double, size_t, size_t, size_t>(), pybind11::arg("numberOfParticles"),
           pybind11::arg("temperature"), pybind11::arg("dt"), pybind11::arg("boxSize"),
           pybind11::arg("sampleFrequency") = 100, pybind11::arg("logLevel") = 0, pybind11::arg("seed") = 12)
      .def_readonly("rdfSampler", &MolecularDynamics::rdfSampler)
      .def_readonly("msdSampler", &MolecularDynamics::msdSampler)
      .def("__repr__", &MolecularDynamics::repr)
      .def("run", &MolecularDynamics::run, pybind11::arg("numberOfSteps"), pybind11::arg("equilibrate"),
           pybind11::arg("outputPDB") = true);
  pybind11::class_<SampleRDF>(m, "SampleRDF").def("getResults", &SampleRDF::getResults);
  pybind11::class_<SampleMSD>(m, "SampleMSD").def("getResults", &SampleMSD::getResults);
}
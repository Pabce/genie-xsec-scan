// Standalone GENIE cross-section scanner.
//
// Unity translation unit for the standalone scanner. Cohesive implementation
// fragments live in detail/ while internal symbols retain one linkage scope.

#include "xsec_scan.hpp"

#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/XclsTag.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/RunOpt.h"
#if __has_include("Framework/Interaction/KPhaseSpaceCuts.h")
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#define XSEC_SCAN_HAS_KPHASESPACECUTS 1
#else
#define XSEC_SCAN_HAS_KPHASESPACECUTS 0
#endif
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace xsec_scan {
namespace {

#include "detail/options.inc"
#include "detail/runtime_points.inc"
#include "detail/kinematics.inc"
#include "detail/qel_folding.inc"
#include "detail/initial_state_folding.inc"
#include "detail/evaluation.inc"
#include "detail/output_run.inc"

} // namespace

int run_cli(int argc, char **argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "xsec_scan: " << e.what() << "\n\n";
    print_usage(std::cerr);
    std::_Exit(1);
  }
}

} // namespace xsec_scan

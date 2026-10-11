// Diagnostics of the global weak (deformation) modes of an alignment result
#ifndef O2ALIGN_WEAKMODES_H
#define O2ALIGN_WEAKMODES_H

#include <map>
#include <string>
#include <cstdint>

namespace o2::alignrs
{
class Volume;

/// Diagnose the global deformation modes of Volume::GlobalMode in the fitted rigid-body corrections.
///
/// For every detector branch (child of the root) two complementary numbers are given per mode:
///  - metric A (what Millepede's mode constraint sees): sum(c_i * p_i) / sum(|c_i|) over the free
///    rigid-body DOFs of the branch, c_i being the coefficient the DOF would have in the mode equation
///    written by Volume::writeGlobalModeConstraints, p_i its fitted value. In micrometers (rotations
///    enter with the lever arm of the mode field), as quoted in AlignmentWeakModes_handover.md.
///  - metric B (what the geometry actually did): the corrections of the whole hierarchy are composed
///    into the displacement of the centre of every sensor, which is fitted by a global rigid shift +
///    rotation plus the mode fields. The amplitudes are given as the displacement the mode produces at
///    the largest lever arm of the branch, with the fraction of the displacement power they explain and
///    the mode-to-mode correlations. Independent of which hierarchy level carries the correction.
/// Writes <outBase>.json (numbers), <outBase>.root (per-sensor displacements and fit results, see
/// macro/drawWeakModes.C) and logs a summary table. The fitted values are absolute: an injected
/// misalignment (misAlgJson) is not subtracted.
void diagnoseWeakModes(Volume* root, const std::map<uint32_t, double>& labelToValue, const std::string& outBase);

} // namespace o2::alignrs

#endif

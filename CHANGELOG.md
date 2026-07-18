# Changelog

## 1.0.2

### Fixed

- Preserved valid decomposition segments when reconciling overlaps between
  long-sequence processing blocks.
- Added regression coverage for exact boundaries, small overlaps, duplicate
  segments, conflicting overlaps, reverse-strand segments, and consecutive
  block seams.

This patch does not change seeding, scoring, wavefront search, backtrace, or
ordinary non-overlap decomposition behavior.

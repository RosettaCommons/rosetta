

#ifndef INCLUDED_protocols_forge_build_SegmentInsert_SegmentInsertConnectionScheme_SegmentInsertConnectionScheme_Enum_hh
#define INCLUDED_protocols_forge_build_SegmentInsert_SegmentInsertConnectionScheme_SegmentInsertConnectionScheme_Enum_hh

namespace protocols {
namespace forge {
namespace build {
// description namespace to hold enum
namespace SegmentInsertConnectionScheme {
	/// @brief connect insertion on its N-side, C-side, or decide randomly
	///  between the two
enum Enum {
	RANDOM_SIDE,
	N,
	C
};

}
}
}
}

#endif

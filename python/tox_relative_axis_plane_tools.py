import numpy as np
import ctypes

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")


def setup_omics_vector_RAP_projection():
    omics_vector_RAP_projection = lib.omics_vector_RAP_projection_c

    omics_vector_RAP_projection.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
    ]

    omics_vector_RAP_projection.restype = None

    def initArgs(vecs, vecs_selection_mask, axes_selection_mask):
        return (
            np.asfortranarray(vecs, dtype=np.float64),
            np.ascontiguousarray(vecs_selection_mask, dtype=np.int32),
            np.ascontiguousarray(axes_selection_mask, dtype=np.int32),
            np.empty((sum(axes_selection_mask), sum(vecs_selection_mask)), order="F", dtype=np.float64)
        )

    return omics_vector_RAP_projection, initArgs


def setup_omics_field_RAP_projection():
    omics_field_RAP_projection = lib.omics_field_RAP_projection_c

    omics_field_RAP_projection.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
    ]

    omics_field_RAP_projection.restype = None

    def initArgs(vecs, vecs_selection_mask, axes_selection_mask):
        return (
            np.asfortranarray(vecs, dtype=np.float64),
            np.ascontiguousarray(vecs_selection_mask, dtype=np.int32),
            np.ascontiguousarray(axes_selection_mask, dtype=np.int32),
            np.empty((sum(axes_selection_mask), sum(vecs_selection_mask)), order="F", dtype=np.float64)
        )

    return omics_field_RAP_projection, initArgs


def test_omics_vector_RAP_projection_call():
    omics_vector_RAP_projection, initArgs = setup_omics_vector_RAP_projection()

    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    vecs = np.random.rand(n_axes, n_vecs)

    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    n_selected_vecs = vecs_selection_mask.sum()

    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    vecs, vecs_selection_mask, axes_selection_mask, projections = initArgs(vecs, vecs_selection_mask, axes_selection_mask)

    # using filled matrix just to make sure it was changed
    projections[:] = 1

    omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)

    for i_vec in range(n_selected_vecs):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"


def test_omics_field_RAP_projection_call():
    omics_field_RAP_projection, initArgs = setup_omics_field_RAP_projection()

    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    vecs = np.random.rand(2 * n_axes, n_vecs)

    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    n_selected_vecs = vecs_selection_mask.sum()

    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    vecs, vecs_selection_mask, axes_selection_mask, projections = initArgs(vecs, vecs_selection_mask, axes_selection_mask)

    # using filled matrix just to make sure it was changed
    projections[:] = 1

    omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)

    for i_vec in range(n_selected_vecs):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"


def main():
    test_omics_vector_RAP_projection_call()
    test_omics_field_RAP_projection_call()
    print("Done")


if __name__ == '__main__':
    main()

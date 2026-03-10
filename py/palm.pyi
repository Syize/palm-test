import numpy as np

class palm:

    terrain_height: np.ndarray

    @classmethod
    def extract_scalar_grid_terrain_height(cls, file_path: str) -> np.ndarray:
        ...

    @classmethod
    def extract_w_grid_terrain_height(cls, file_path: str) -> np.ndarray:
        ...

    @classmethod
    def say_hello(cls) -> None:
        ...
    
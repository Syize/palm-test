"""
Interface to read terrain height data in PALM.
"""

import numpy as np

# extension should be place in the same directory as this script.
from . import palm_model


class PALMTopoExtracter:
    __initialized = False
    __instance = None
    
    def __init__(self, file_path: str) -> None:
        if len(file_path) > 100:
            raise ValueError("`file_path` can't be longer than 100.")
        
        self.file_path = file_path

    def __new__(cls, *args) -> "PALMTopoExtracter":
        if not cls.__initialized or cls.__instance is None:
            cls.__initialized = True
            cls.__instance = super().__new__(cls)

        else:
            print("If you want to read new static driver, re-run your python script.")
            print("INIT EXTRACTER TWICE WON'T WORK!!!")

        return cls.__instance

    def read_scalar_terrain_height(self) -> np.ndarray:
        palm_model.palm_test.extract_scalar_grid_terrain_height(self.file_path)
    
        terrain_height = palm_model.palm_test.terrain_height.copy()
    
        return terrain_height
    
    def read_w_terrain_height(self) -> np.ndarray:
        palm_model.palm_test.extract_w_grid_terrain_height(self.file_path)
    
        terrain_height = palm_model.palm_test.terrain_height.copy()
    
        return terrain_height


__all__ = ["PALMTopoExtracter"]

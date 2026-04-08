import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

def ij2n (i, j, N):
    return i + j*N


{-# LANGUAGE MultiParamTypeClasses, AllowAmbiguousTypes, GADTs #-}

import Data.Map (Map)
import qualified Data.Map as Map

class MarkovCategory o m where
    dom :: m -> o
    cod :: m -> o
    sharp :: m -> m
    samp :: o -> m
    cond :: m -> [o] -> m
    id :: o -> m
    arrange :: [Int] -> [m] -> Maybe m
    compose :: m -> m -> Maybe m

data FinStoch p where
    distr :: Map a p -> FinStoch p

import pywt
import numpy as np

class WaveletEMD:
    """
    Wavelet EMD using multilevel wavelet transform of 1d-vector
    """

    def __init__(self, wavelet: str = "sym5", level: int = 6, mode: str = "zero"):
        self.level = level
        self.mode  = mode 
        self.wavelet = wavelet

    def __call__(self, vec1: np.ndarray, vec2: np.ndarray) -> float:
        """Compute waveletEMD between two histogram descriptors

        Args:
            vec1 (np.ndarray): histogram descriptor 1
            vec2 (np.ndarray): histogram descriptor 2

        Returns:
            float: waveletEMD between vec1 and vec2
        """        
        # compute wavelet coefficients for each vector
        coeffs1, coeffs2 = self.to_wavelet_domain(vec1), self.to_wavelet_domain(vec2)
        
        # rescale coefficients (ommit A1)
        res_coeffs1, res_coeffs2 = self.rescale_coefficients(coeffs1), self.rescale_coefficients(coeffs2)

        # compute distance
        dist = np.abs(np.concatenate(res_coeffs1) - np.concatenate(res_coeffs2)).sum()

        return dist


    def to_wavelet_domain(self, vec: np.ndarray) -> np.ndarray:
        "Compute wavelet coefficients for a histogram descriptor"
        # compute wavelet coefficients  
        coeffs = pywt.wavedec(data=vec, wavelet=self.wavelet, mode=self.mode, level=self.level)

        return coeffs

    def rescale_coefficients(self, coeffs: np.ndarray) -> np.ndarray:
        """ 
        Rescale coefficients as shown in 
        Eq (2) doi: 10.1109/CVPR.2008.4587662
        """
                
        j = self.level # starting level
        s = 1 # l1 (Manhattan) distance
        n = 1 # dimension

        # rescale coefficients
        res_coeffs = []
        for coeff in coeffs[1:]:
            scale_factor = 2**(-j*(s+n/2))
            res_coeffs.append(
                scale_factor * coeff
            )
            j -= 1 # update level
        return res_coeffs
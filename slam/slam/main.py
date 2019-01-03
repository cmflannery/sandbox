import os
import time
import numpy as np

import cv2
from skimage import feature
import pygame
import fire

debug = True

W = 1280//1
H = 640//1

orb = cv2.ORB_create()

fast = cv2.FastFeatureDetector_create()

def slam(frame, disp):
    """Executes the SLAM algorithm.

    Parameters
    ----------
    frame : np.ndarray
        current frame, expected rgb

    Returns
    -------
    tbd
    """
    kp1 = cv2.goodFeaturesToTrack(np.mean(frame,axis=2).astype(np.uint8), 3000, qualityLevel=0.01, minDistance=3)
    for p in kp1:
        u,v = map(lambda x: int(round(x)), p[0])
        cv2.circle(frame, (u,v), color=(0,255,0), radius=3)

    print(frame.shape)
    disp.paint(frame)


class FeatureExtractor(object):
    GX = 16
    GY = 16

    def __init__(self):
        pass

    def extract(self, image):
        return cv2.goodFeaturesToTrack(image,qualityLevel=0.01)


class Display():
    def __init__(self, W, H):
        pygame.init()
        self.screen = pygame.display.set_mode((W,H))


    def paint(self,image):
        surf = pygame.surfarray.make_surface(image.swapaxes(0,1)).convert()
        self.screen.blit(surf, (0,0))
        pygame.display.flip()
        pygame.display.update()

fe = FeatureExtractor()

def main(file_name):
    file_path = os.path.join(os.path.dirname(__file__),file_name)

    cap = cv2.VideoCapture(file_path)
    disp = Display(W,H)
    overflow = 0
    while cap.isOpened() and overflow < 100:
        overflow += 0
        retval, frame = cap.read()

        if retval:
            slam(frame, disp)
        else:
            continue


    print(cap)

if __name__ == "__main__":
    fire.Fire(main)

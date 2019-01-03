#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


def get_data():
    with open('3le.txt', 'r') as f:
        data = f.read()
    return data

def parse_data(data):

    satellites = {}
    for line in data.splitlines():
        print(line[0])


def main():
    parse_data(get_data())


if __name__ == "__main__":
    main()

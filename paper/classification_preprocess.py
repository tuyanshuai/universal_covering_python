import numpy as np
import os

categories = ['bathtub', 'bed', 'chair', 'desk', 'dresser', 'monitor', 'night_stand', 'sofa', 'table', 'toilet']
path = '/home/local/ASUAD/yanshuai/shapeAI/Data/yanshuai/ModelNet10/'
pclexe = 'pcl_mesh_sampling'

def OFFtoPLY(path, categories, DataGroup):
    for cat in categories:
        DataArray = []
        # deal with train first
        files = os.listdir(path + cat + '/' + DataGroup + '/')
        files = [x for x in files if x[-4:] == '.off']
        for file_index, file in enumerate(files):
            fileName = file.split('.')[0]
            with open(path + cat + '/' + DataGroup + '/' + file, 'r') as f:
                tmp = f.readline().replace('\n', '')
                line = ''
                if tmp != 'OFF':
                    line = tmp[3:]
                else:
                    line = f.readline().replace('\n', '')

                # get number of points in the model
                point_count = line.split(' ')[0]
                face_count = line.split(' ')[1]

                data = []
                # fill ndarray with datapoints
                for index in range(0, int(point_count)):
                    line = f.readline().rstrip().split()
                    line[0] = float(line[0])
                    line[1] = float(line[1])
                    line[2] = float(line[2])
                    data.append(line)
                data = np.array(data)
                # normalize data before conversion
                centroid = np.mean(data, axis=0)
                data = data - centroid
                m = np.max(np.sqrt(np.sum(data ** 2, axis=1)))
                data = data / m

                # create ply file,write in header first.
                with open(path + cat + '/' + DataGroup + '/' + fileName + ".ply", 'w') as plyFile:
                    plyFile.write('ply\nformat ascii 1.0\nelement vertex ')
                    plyFile.write(point_count)
                    plyFile.write('\nproperty float32 x\nproperty float32 y\nproperty float32 z\nelement face ')
                    plyFile.write(face_count)
                    plyFile.write('\nproperty list uint8 int32 vertex_indices\nend_header\n')
                    for index in range(0, int(point_count)):
                        plyFile.write(' '.join(map(str, data[index])))
                        plyFile.write('\n')
                    for index in range(0, int(face_count)):
                        plyFile.write(f.readline())


import subprocess


def PLYtoPCD(path, categories, DataGroup):
    for cat in categories:
        DataArray = []
        # deal with train first
        files = os.listdir(path + cat + '/' + DataGroup + '/')
        files = [x for x in files if x[-4:] == '.ply']
        for file_index, file in enumerate(files):
            fileName = file.split('.')[0]
            subprocess.call([pclexe, path + cat + '/' + DataGroup + '/' + file,
                             path + cat + '/' + DataGroup + '/' + fileName + ".pcd", '-no_vis_result', '-n_samples',
                             '200000', '-leaf_size', '0.01'])


def delPCDHeader(path, categories, DataGroup):
    for cat in categories:
        files = os.listdir(path + cat + '/' + DataGroup + '/')
        files = [x for x in files if x[-4:] == '.pcd']
        for file_index, file in enumerate(files):
            with open(path + cat + '/' + DataGroup + '/' + file, 'r') as f:
                for y in range(9):
                    f.readline()
            with open(path + cat + '/' + DataGroup + '/' + file, 'r') as f:
                data = f.read().splitlines(True)
            with open(path + cat + '/' + DataGroup + '/' + file, 'w') as f:
                f.writelines(data[9:])


if __name__ == "__main__":
    # OFFtoPLY(path, categories, 'train')
    # OFFtoPLY(path, categories, 'test')
    #
    # PLYtoPCD(path, categories, 'train')
    # PLYtoPCD(path, categories, 'test')

    delPCDHeader(path, categories, 'test')
    delPCDHeader(path, categories, 'test')
clear; close all; clc;

%% 2次元ベクトルの可視化
% 2次元空間上のベクトルを定義
u1 = [2, 1]';
u2 = [-3, 2]';
u3 = [-1, -3]';
u4 = 2*u1;
u5 = -0.5*u3;

% 行列Uとしてまとめる（各列がベクトル）
U = [u1, u2, u3, u4, u5];

% ベクトルを可視化
plot2d(U);

%% 2次元空間上の部分空間の可視化
% u2が生成する部分空間を可視化
plot2d(U, 2);

% u1及びu4が生成する部分空間(=u1が生成する部分空間)を可視化
plot2d(U, [1, 4]);

% u1及びu2が生成する部分空間(=2次元空間)を可視化
plot2d(U, [1, 2]);

%% 3次元ベクトルの可視化
% 3次元空間上の3個のベクトルを定義
u1 = [1, -3, 0]';
u2 = [2, 0, -2]';
u3 = [3, 3, 3]';
u4 = 2*u1;
u5 = [-2, 0, 5]';

% 行列Uとしてまとめる（各列がベクトル）
U = [u1, u2, u3, u4, u5];

% ベクトルを可視化
plot3d(U);

%% 3次元空間上の部分空間の可視化
% u3が生成する部分空間を可視化
plot3d(U, 3);

% u1及びu5が生成する部分空間を可視化
plot3d(U, [1, 5]);

% u1及びu4が生成する部分空間(=u1が生成する部分空間)を可視化
plot3d(U, [1, 4]);

% u2及びu5が生成する部分空間(=xz平面)を可視化
plot3d(U, [2, 5]);

% u1, u2, 及びu3が生成する部分空間(=3次元空間)を可視化
plot3d(U, [1, 2, 3]);
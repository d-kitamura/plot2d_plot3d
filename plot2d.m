function plot2d(vec, subSp)
% plot2d: 入力ベクトルを2次元ベクトル空間上に矢印として描画
%
% [Input]
%   vecs: 入力ベクトルを各列に持つ行列 ([u1], [u1, u2], [u1, u2, u3], ...)
%  subSp: 1個又は複数個のベクトルが生成する部分空間の描画
%         スカラーでnを与えるとn個目のベクトルの部分空間を描画
%         配列で[n1, n2]を与えるとn1個目とn2個目のベクトルが生成する部分空間を描画
% 
arguments
    vec (2,:) {mustBeNumeric}
    subSp (1, :) {mustBeInteger} = []
end

nVec = size(vec, 2); % 入力されたベクトルの本数
lw = 0.75; % 線の太さ
m = max(abs(vec), [], "all"); % 原点から最も離れている座標

% 2次元空間にベクトルを描画
figure; hold on;
for iVec = 1:nVec
    quiver(0, 0, vec(1, iVec), vec(2, iVec), "AutoScale", "off", "LineWidth", lw+1, "MaxHeadSize", 0.3*m/norm(vec(:, iVec))); % 矢印の描画
    lgd(iVec) = "(" + vec(1, iVec) + ", " + vec(2, iVec) + ")"; % 凡例
end
ax = gca;
box on; grid on; ax.LineWidth = lw;
minEnd = -m-2; maxEnd = m+2; % 軸の表示範囲の最小値と最大値
xlim([minEnd, maxEnd]); ylim([minEnd, maxEnd]);
xticks(minEnd*2:maxEnd*2); yticks(minEnd*2:maxEnd*2);
xlabel("x"); ylabel("y");
plot3(ax.XLim, [0, 0], [0, 0], "k", "LineWidth", lw); % x軸の直線
plot3([0, 0], ax.YLim, [0, 0], "k", "LineWidth", lw); % y軸の直線
legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示

% 平行なベクトルは短いものを前面，長いものを背面に入れ替え
parIdx = pickParallel(vec); % 平行なベクトルの確認
for iIdx = 1:size(parIdx, 1)
    i1 = parIdx(iIdx, 1); i2 = parIdx(iIdx, 2);
    parVec = vec(:, [i1, i2]); % 平行なベクトル2個をvecから抽出
    if norm(parVec(:, 2)) > norm(parVec(:, 1)) % 後から描かれた（前面にある）ベクトルの方が長いとき
        hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得
        tmp = hdl(end-(i1-1)); % 短いベクトルと
        hdl(end-(i1-1)) = hdl(end-(i2-1)); % 長いベクトルの
        hdl(end-(i2-1)) = tmp; % 前面と背面を入れ替え
        set(ax, "children", hdl); % 軸の子クラスを変更後で更新
    end
end

% 部分空間の描画
if exist("subSp", "var") && ~isempty(subSp) % subSpが引数として与えられているとき
    r = rank(vec(:, subSp)); % 指定されたベクトルの1次独立な最大個数を取得
    if r == 1 % 部分空間は直線（1次元空間）
        u = vec(:, subSp(1));
        k = minEnd:maxEnd; % 媒介変数（k*uの係数k）
        X = k*u(1); Y = k*u(2);
        plot(X, Y, "Color", [0.85, 0.85, 0.85], "LineWidth", lw+6);
        legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示（この後順番を入れ替えるのでここで表示しておく）
        hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得
        hdl = hdl([2:end, 1]); % 最新の描画（部分空間のplot3）を一番背面に変更
        set(ax, "children", hdl); % 軸の子クラスを変更後で更新
    else % 部分空間は平面（2次元空間）
        set(ax, "Color", [0.85, 0.85, 0.85]); % figureの背景を透明に設定
    end
else
    legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示
end

% x軸とy軸の直線を最背面に移動
hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得
hdl = hdl([3:end, 1, 2]); % 前面・背面を入れ替え
set(ax, "children", hdl); % 軸の子クラスを変更後で更新
end

% -------------------------------------------------------------------------
function idx = pickParallel(vec)
ptn = nchoosek(1:size(vec, 2), 2); % n個のベクトルから2個選ぶ組み合わせの全パターン
nPtn = size(ptn, 1); % パターンの総数
r = zeros(nPtn, 1);
for iPtn = 1:nPtn
    idx = ptn(iPtn, :); % チェックする2個のベクトルのインデクス
    r(iPtn) = rank(vec(:, idx)); % 2個のベクトルを列に持つ行列のランクを調べる（1なら1次従属）
end
idx = ptn(r<2, :); % 1次従属なベクトルの組み合わせ
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
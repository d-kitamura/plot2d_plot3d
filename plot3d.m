function plot3d(vec, subSp)
% plot3d: 入力ベクトルを3次元ベクトル空間上に矢印として描画
%         subSpでベクトルを指定することで部分空間の描画も可能
%
% [Input]
%    vec: 入力ベクトルを各列に持つ行列 ([u1], [u1, u2], [u1, u2, u3], ...)
%  subSp: 1個又は複数個のベクトルが生成する部分空間の描画
%         スカラーでnを与えるとn個目のベクトルの部分空間を描画
%         配列で[n1, n2]を与えるとn1個目とn2個目のベクトルが生成する部分空間を描画
%
arguments
    vec (3, :) {mustBeNumeric}
    subSp (1, :) {mustBeInteger, mustBePositive} = []
end

nVec = size(vec, 2); % 入力されたベクトルの本数
lw = 0.75; % 線の太さ
m = max(abs(vec), [], "all"); % 原点から最も離れている座標

% 3次元空間にベクトルを描画
figure; hold on;
for iVec = 1:nVec
    quiver3(0, 0, 0, vec(1, iVec), vec(2, iVec), vec(3, iVec), "AutoScale", "off", "LineWidth", lw+1, "MaxHeadSize", 0.4*m/norm(vec(:, iVec))); % 矢印の描画
    lgd(iVec) = "(" + vec(1, iVec) + ", " + vec(2, iVec) + ", " + vec(3, iVec) + ")"; % 凡例
end
ax = gca;
box on; grid on; ax.LineWidth = lw;
ax.YDir = 'reverse'; % y軸は手前が正
minEnd = -m-2; maxEnd = m+2; % 軸の表示範囲の最小値と最大値
xlim([minEnd, maxEnd]); ylim([minEnd, maxEnd]); zlim([minEnd, maxEnd]);
xticks(minEnd*2:maxEnd*2); yticks(minEnd*2:maxEnd*2); zticks(minEnd*2:maxEnd*2);
plot3(ax.XLim, [0, 0], [0, 0], "k", "LineWidth", lw); % x軸の直線
plot3([0, 0], ax.YLim, [0, 0], "k", "LineWidth", lw); % y軸の直線
plot3([0, 0], [0, 0], ax.ZLim, "k", "LineWidth", lw); % z軸の直線
xlabel("x"); ylabel("y"); zlabel("z");

% 平行なベクトルは短いものを前面，長いものを背面に入れ替え
parIdx = pickParallel(vec); % 平行なベクトルの確認
for iIdx = 1:size(parIdx, 1)
    lgd = swapVecs(parIdx(iIdx, :), vec, ax, lgd);
end

% 部分空間の描画
if exist("subSp", "var") && ~isempty(subSp) % subSpが引数として与えられているとき
    subSp = subSp(sum(vec(:, subSp)~=0, 1)~=0); % 零ベクトルは部分空間に寄与しないので除外
    r = rank(vec(:, subSp)); % 指定されたベクトルの1次独立な最大個数を取得
    if r == 0 % 部分空間は原点（0次元空間）
        scatter3(0, 0, 0, 36, "filled", "MarkerFaceColor", [0.85, 0.85, 0.85], "MarkerEdgeColor", [0.85, 0.85, 0.85]);
        legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示（この後順番を入れ替えるのでここで表示しておく）
    elseif r == 1 % 部分空間は直線（1次元空間）
        u = vec(:, subSp(1));
        k = minEnd:maxEnd; % 媒介変数（k*uの係数k）
        X = k*u(1); Y = k*u(2); Z = k*u(3);
        plot3(X, Y, Z, "Color", [0.85, 0.85, 0.85], "LineWidth", lw+6);
        legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示（この後順番を入れ替えるのでここで表示しておく）
        hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得
        hdl = hdl([2:end, 1]); % 最新の描画（部分空間のplot3）を一番背面に変更
        set(ax, "children", hdl); % 軸の子クラスを変更後で更新
    elseif r == 2 % 部分空間は平面（2次元空間）
        u = vec(:, subSp(1));
        v = vec(:, subSp(2));
        w = cross(u, v); % uとvが生成する平面の法線ベクトル（外積）
        [X, Y] = meshgrid(minEnd:maxEnd, minEnd:maxEnd); % xy平面のメッシュグリッド定義
        if w(3) == 0; w(3) = 10^-10; end % 部分空間がxy平面に対して垂直になるときの0割りを回避
        Z = -w(1)/w(3)*X + -w(2)/w(3)*Y; % 法線ベクトルを使った平面方程式
        sObj = surf(X, Y, Z, "FaceAlpha", 0.5, "FaceColor", [0.7, 0.7, 0.7]);
        sObj.EdgeColor = 'none';
        legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示（この後順番を入れ替えるのでここで表示しておく）
        hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得
        hdl = hdl([2:end, 1]); % 最新の描画（部分空間のsurf）を一番背面に変更
        set(ax, "children", hdl); % 軸の子クラスを変更後で更新
    else % 部分空間は3次元空間
        set(ax, "Color", [0.85, 0.85, 0.85]); % figureの背景を灰色に設定
        set(ax, "BoxStyle", "full");
        legend(lgd, "Location", "southwest", "LineWidth", lw, "Color", "w"); % 凡例の表示（背景色は再び白に変更）
    end
else
    legend(lgd, "Location", "southwest", "LineWidth", lw); % 凡例の表示
end
view(40, 20);
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

% -------------------------------------------------------------------------
function lgd = swapVecs(idx, vec, ax, lgd)
    i1 = idx(1); i2 = idx(2);
    parVec = vec(:, [i1, i2]); % 平行なベクトル2個をvecから抽出
    if norm(parVec(:, 2)) > norm(parVec(:, 1)) && norm(parVec(:, 1)) ~= 0 % 後から描かれた（前面にある）ベクトルの方が長いとき
        hdl = get(ax, "Children"); % 軸の子クラスから描画の順番を取得

        % 短いベクトルと長いベクトルの前面と背面を入れ替え
        tmp = hdl(end-(i1-1));
        hdl(end-(i1-1)) = hdl(end-(i2-1));
        hdl(end-(i2-1)) = tmp;
        set(ax, "children", hdl); % 軸の子クラスを変更後で更新

        % 短いベクトルと長いベクトルの凡例も入れ替え
        tmp = lgd(i1);
        lgd(i1) = lgd(i2);
        lgd(i2) = tmp;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
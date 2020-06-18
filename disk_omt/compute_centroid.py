"""
cc = compute_centroid(cp,pd)


    cc = zeros(length(pd.cell), 2);
    for i = 1:length(pd.cell)
    ci = pd.cell
    {i};
    di = pd.dpe(ci,:);
    if (isempty(ci))
        warning('empty ci');
    cc(i,:) = [0 0];
    else

    if ci(1) == ci(end)
        di = di(1:end - 1,:);
        end

        pc = polybool([cp(:, 1), cp(:, 2)], [di(:, 1), di(:, 2)], 'and', [], [], 10 ^ 10);
        if length(pc) ~= 1
        % pause
        warning('length(pc) ~= 1');

        cc(i,:)= mean(di) / norm(mean(di));
        continue;
    end

    cc(i,:) = centroid(pc
    {1});
    end


end
"""

import numpy as np


def compute_centroid(cp, pd):
    cc = np.zeros( (pd["cell"].shape[0],2 ) )


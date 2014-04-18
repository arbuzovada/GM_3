function energy = get_energy_ab(alpha, beta, labels, unary_pot, intensity)
% This function calculates energy for given label configuration
% INPUT:
%    labels: N-by-M array of integer, labels for pixels
%    unary_pot: N-by-M-by-K array of double, unary potentials
%    pair_pot: N-by-M-by-K array of double, pixel brightness
%
% OUTPUT:
%    energy: double

    [N, M] = size(labels);
    energy = 0;
    for i = 1 : N
        for j = 1 : M
            if labels(i, j) == alpha || labels(i, j) == beta
                energy = energy + unary_pot(i, j, labels(i, j));
                if (i < N)
                    energy = energy + get_pp(i + N * (j - 1), i + 1 + N * (j - 1), ...
                        labels(i, j), labels(i + 1, j), intensity);
                end
                if (j < M)
                    energy = energy + get_pp(i + N * (j - 1), i + N * j, ...
                        labels(i, j), labels(i, j + 1), intensity);
                end
            else
                if (i < N)
                    if labels(i + 1, j) == alpha || labels(i + 1, j) == beta
                        energy = energy + get_pp(i + N * (j - 1), i + 1 + N * (j - 1), ...
                         labels(i, j), labels(i + 1, j), intensity);
                    end
                end
                if (j < M)
                    if labels(i, j + 1) == alpha || labels(i, j + 1) == beta
                        energy = energy + get_pp(i + N * (j - 1), i + N * j, ...
                        labels(i, j), labels(i, j + 1), intensity);
                    end
                end
            end
        end
    end
end
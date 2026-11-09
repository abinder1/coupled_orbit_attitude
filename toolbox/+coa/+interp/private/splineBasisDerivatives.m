% Copyright (c) 2026 Andrew Binder

function ders = splineBasisDerivatives(i, u, p, n, knots)
    % SPLINEBASISDERIVATIVES Nonzero B-spline basis derivatives at u.
    %   i is the zero-based knot span, p is the degree, and n is the
    %   maximum derivative order. Knots are stored in a one-based array.
    %   Algorithm A2.3 from Piegl and Tiller, The NURBS Book.

    ndu = zeros(p + 1);
    left = zeros(1, p + 1);
    right = zeros(1, p + 1);
    ndu(1, 1) = 1;

    for j = 1:p
        left(j + 1) = u - knots(i + 2 - j);
        right(j + 1) = knots(i + j + 1) - u;
        saved = 0;
        for r = 0:j - 1
            ndu(j + 1, r + 1) = right(r + 2) + left(j - r + 1);
            temp = ndu(r + 1, j) / ndu(j + 1, r + 1);
            ndu(r + 1, j + 1) = saved + right(r + 2) * temp;
            saved = left(j - r + 1) * temp;
        end
        ndu(j + 1, j + 1) = saved;
    end

    ders = zeros(n + 1, p + 1);
    ders(1, :) = ndu(:, p + 1).';
    a = zeros(2, p + 1);
    for r = 0:p
        s1 = 1;
        s2 = 2;
        a(1, 1) = 1;
        for kk = 1:n
            dv = 0;
            rk = r - kk;
            pk = p - kk;
            if r >= kk
                a(s2, 1) = a(s1, 1) / ndu(pk + 2, rk + 1);
                dv = a(s2, 1) * ndu(rk + 1, pk + 1);
            end
            if rk >= -1
                j1 = 1;
            else
                j1 = -rk;
            end
            if r - 1 <= pk
                j2 = kk - 1;
            else
                j2 = p - r;
            end
            for j = j1:j2
                a(s2, j + 1) = (a(s1, j + 1) - a(s1, j)) / ndu(pk + 2, rk + j + 1);
                dv = dv + a(s2, j + 1) * ndu(rk + j + 1, pk + 1);
            end
            if r <= pk
                a(s2, kk + 1) = -a(s1, kk) / ndu(pk + 2, r + 1);
                dv = dv + a(s2, kk + 1) * ndu(r + 1, pk + 1);
            end
            ders(kk + 1, r + 1) = dv;
            [s1, s2] = deal(s2, s1);
        end
    end

    factor = p;
    for kk = 1:n
        ders(kk + 1, :) = ders(kk + 1, :) * factor;
        factor = factor * (p - kk);
    end
end

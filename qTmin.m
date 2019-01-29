function [qT, ts] = qTmin(Y,Yref, delta, i_R_peaks, iBefore, iAfter,p, size_N,K_starter) 

if K_starter > 2
    p = p + 1;
end
Q = [];
min = NaN;
mindelta = NaN;
epsilon = NaN;
Best_Q = [];

jta = zeros(iAfter+iBefore+1);

    for i=-delta:delta
        upperzero = zeros((delta + i), iAfter+iBefore+1);
        lowerzero = zeros ((delta - i), iAfter+iBefore+1);
        IMatrix = eye (iAfter+iBefore+1);
        jta_t = [upperzero' IMatrix' lowerzero'];
        jta = jta_t';
        
        if i_R_peaks(p) + iAfter > size(Y,2)
            zeta = Y(:, i_R_peaks(p) - iBefore : size(Y,2)) *jta' * Yref';
        else
            zeta = Y(:, i_R_peaks(p) - iBefore : i_R_peaks(p) + iAfter ) *jta' * Yref';
        end
        
        [U,S,V] = svd (zeta);
        Q = U * V';
        if i_R_peaks(p) + iAfter > size(Y,2)
            epsilon = norm(Y(:, i_R_peaks(p) - iBefore : size(Y,2)) - Q * Yref * jta, 'fro');
        else
            epsilon = norm(Y(:, i_R_peaks(p) - iBefore : i_R_peaks(p) + iAfter ) - Q * Yref * jta, 'fro');
        end
        if i == -delta %initial first min and mindelta
            min = epsilon;
            mindelta = i;
            Best_Q = Q;
 
        end
        if epsilon < min
            min = epsilon;
            mindelta = i;
            Best_Q = Q;
            
        end
    end
    qT = Best_Q;
    ts = mindelta;
    
plot(Y(1, i_R_peaks(p) - iBefore: i_R_peaks(p) + iAfter))
hold on
plot(Yref(1,:));
hold on
ym = Best_Q * Yref * jta;
plot(ym(1,:),'m');
hold off


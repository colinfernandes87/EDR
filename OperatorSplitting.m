function [optimal_Q, Ts]=OperatorSplitting(Y,Yref, delta, i_R_peaks, iBefore, iAfter,p, size_N, K_starter)

if K_starter > 2
    p = p + 1;
end
first_Q = eye(3);
min = NaN;
mindelta = NaN;
epsilon = NaN;
jta = zeros(iAfter+iBefore+1);
for k = 1:3
    for i=-delta:delta
        upperzero = zeros((delta + i), iAfter+iBefore+1);
        lowerzero = zeros ((delta - i), iAfter+iBefore+1);
        IMatrix = eye (iAfter+iBefore+1);
        jta_t = [upperzero' IMatrix' lowerzero'];
        jta = jta_t';
        if i_R_peaks(p) + iAfter > size(Y,2)
            epsilon = norm(Y(:, i_R_peaks(p) - iBefore : size(Y,2)) - first_Q * Yref * jta, 'fro');
        else
            epsilon = norm(Y(:, i_R_peaks(p) - iBefore : i_R_peaks(p) + iAfter ) - first_Q * Yref * jta, 'fro');
        end
        
        if i == -delta %initial first min and mindelta
            min = epsilon;
            mindelta = i;
            
        end
        if epsilon < min
            min = epsilon;
            mindelta = i;
            
        end
    end

    upperzero = zeros((delta + mindelta), iAfter+iBefore+1);
    lowerzero = zeros ((delta - mindelta), iAfter+iBefore+1);
    IMatrix = eye (iAfter+iBefore+1);
    jta_t = [upperzero' IMatrix' lowerzero'];
    jta = jta_t';
    if i_R_peaks(p) + iAfter > size(Y,2)
        zeta = Y(:, i_R_peaks(p) - iBefore : size(Y,2)) *jta' * Yref';
    else 
        zeta = Y(:, i_R_peaks(p) - iBefore : i_R_peaks(p) + iAfter ) *jta' * Yref';
    end
    
    [U,S,V] = svd (zeta);
    first_Q = U * V';
    
end
optimal_Q = first_Q;
Ts = mindelta;

plot(Y(1, i_R_peaks(p) - iBefore: i_R_peaks(p) + iAfter))
hold on
plot(Yref(1,:));
hold on
ym = optimal_Q * Yref * jta;
plot(ym(1,:),'m');
hold off

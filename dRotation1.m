function C = dRotation1(n, mu, m, theta)
%%生成旋转系数
% theta = theta * pi / 180; %转化为弧度制\
% 
% log_factor = 0.5 * (gammaln(abs(n + mu) + 1) + gammaln(abs(n - mu) + 1) ...
%                   - gammaln(abs(n + m) + 1) - gammaln(abs(n - m) + 1));
% p1 = exp(log_factor);
% 
% p2 = (cos(theta/2) .^ (mu + m)) * (sin(theta/2) .^ (mu + m));
% C = p1 * p2 * jacobiP(abs(n - mu), abs(mu - m), abs(mu + m), cos(theta));
% end
% 计算阶乘项的系数（使用gammaln避免溢出）
    if theta == 0
        C = (mu == m);
        return
    elseif theta == pi
        C = (-1) ^ (n + m) * (mu == -m);
        return
    end
    log_factor = 0.5 * (gammaln(n + mu + 1) + gammaln(n - mu + 1) ...
                      - gammaln(n + m + 1) - gammaln(n - m + 1));
    factor = exp(log_factor);
    
    % 计算三角函数项
    cos_theta_half = cos(theta / 2);
    sin_theta_half = sin(theta / 2);
    cos_term = cos_theta_half ^ (mu + m);
    sin_term = sin_theta_half ^ (mu - m);
    
    % 计算雅可比多项式
    k = n - mu;       % 多项式次数
    a = mu - m;       % 参数alpha
    b = mu + m;       % 参数beta
    x = cos(theta);   % 变量
    P = jacobiP(k, a, b, x);
    
    % 组合所有项
    D = factor * cos_term * sin_term * P;
    C = D;
end
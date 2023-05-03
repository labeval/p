
clc
clear all

%         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |
linedata = [1      2       0.01938   0.05917    0.0264         1
             1      5       0.05403   0.22304    0.0246         1
             2      3       0.04699   0.19797    0.0219         1
             2      4       0.05811   0.17632    0.0170         1
             2      5       0.05695   0.17388    0.0173         1
             3      4       0.06701   0.17103    0.0064         1
             4      5       0.01335   0.04211    0.0            1
             4      7       0.0       0.20912    0.0        0.978
             4      9       0.0       0.55618    0.0        0.969
             5      6       0.0       0.25202    0.0        0.932
             6     11       0.09498   0.19890    0.0            1
             6     12       0.12291   0.25581    0.0            1
             6     13       0.06615   0.13027    0.0            1
             7      8       0.0       0.17615    0.0            1
             7      9       0.0       0.11001    0.0            1
             9     10       0.03181   0.08450    0.0            1
             9     14       0.12711   0.27038    0.0            1
            10     11       0.08205   0.19207    0.0            1
            12     13       0.22092   0.19988    0.0            1
            13     14       0.17093   0.34802    0.0            1 ];

%         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi |  Qmin | Qmax |
busd =     [1     1    1.060   0       0     0     0     0       0       0;
            2     2    1.045   0      40   42.4  21.7   12.7    -40     50;
            3     2    1.010   0       0   23.4  94.2   19.0     0      40;
            4     3    1.0     0       0     0   47.8   -3.9     0       0;
            5     3    1.0     0       0     0    7.6    1.6     0       0;
            6     2    1.070   0       0   12.2  11.2    7.5    -6      24;
            7     3    1.0     0       0     0    0.0    0.0     0       0;
            8     2    1.090   0       0   17.4   0.0    0.0    -6      24;
            9     3    1.0     0       0     0   29.5   16.6     0       0;
            10    3    1.0     0       0     0    9.0    5.8     0       0;
            11    3    1.0     0       0     0    3.5    1.8     0       0;
            12    3    1.0     0       0     0    6.1    1.6     0       0;
            13    3    1.0     0       0     0   13.5    5.8     0       0;
            14    3    1.0     0       0     0   14.9    5.0     0       0;];
  
fb = linedata(:,1);     % From bus number...
tb = linedata(:,2);     % To bus number...
r = linedata(:,3);      % Resistance, R...
x = linedata(:,4);      % Reactance, X...
b = linedata(:,5);      % Ground Admittance, B/2...
a = linedata(:,6);      % Tap setting value..
z = r + i*x;            % Z matrix...
y = 1./z;               % To get inverse of each element...
b = i*b;                % Make B imaginary...
nb = max(max(fb),max(tb));    % no. of buses...
nl = length(fb);           % no. of branches...
Y = zeros(nb,nb);        % Initialise YBus...
 
 % Formation of the Off Diagonal Elements...
 for k = 1:nl
     Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end
 
 % Formation of Diagonal Elements....
 for m = 1:nb
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end
 Y  ;                % Bus Admittance Matrix
 %zbus = inv(Y);      % Bus Impedance Matrix




%num = 14;   % IEEE-14, IEEE-30, IEEE-57, 118..
%Y = ybusppg(num);           % Calling ybusppg.m to get Bus Admittance Matrix..
%busd = busdatas(num);       % Calling busdata30.m to get busdatas..
BMva = 100;                 % Base MVA..
bus = busd(:,1) ;           % Bus Number..
type = busd(:,2);           % Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,3);              % Specified Voltage..
del = zeros(length(V),1);   % Voltage Angle..
Pg = busd(:,5)/BMva;        % PGi..
Qg = busd(:,6)/BMva;        % QGi..
Pl = busd(:,7)/BMva;        % PLi..
Ql = busd(:,8)/BMva;        % QLi..
Qmin = busd(:,9)/BMva;      % Minimum Reactive Power Limit..
Qmax = busd(:,10)/BMva;     % Maximum Reactive Power Limit..
nbus = max(bus);            % To get no. of buses..
P = Pg - Pl;                % Pi = PGi - PLi..
Q = Qg - Ql;                % Qi = QGi - QLi..
Psp = P;                    % P Specified       
Qsp = Q;                    % Q Specified
G = real(Y);                % Conductance..
B = imag(Y);                % Susceptance..
pv = find(type == 2 | type == 1);       % Index of PV Buses..
pq = find(type == 3);                   % Index of PQ Buses..
npv = length(pv)                       % Number of PV buses..
npq = length(pq);                       % Number of PQ buses..
Er = 1;  
Iter = 1;           % iteration starting
while (Er > 1e-5)   % Iteration starting..
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculate P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
    % Checking Q-limit violations..
    if Iter <= 7 && Iter > 2    % Only checked up to 7th iterations..
        for n = 2:nbus
            if type(n) == 2
                QG = Q(n)+Ql(n);
                if QG < Qmin(n)
                    V(n) = V(n) + 0.01;
                elseif QG > Qmax(n)
                    V(n) = V(n) - 0.01;
                end
            end
         end
    end
    
    % Calculate change from specified value
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    J4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    J = [J1 J2; J3 J4];         % Jacobian Matrix
    
    X = J\M;                    % INV(J) x M, Correction Vector..
    dTh = X(1:nbus-1);          % Change in Voltage Angle..
    dV = X(nbus:end);           % Change in Voltage Magnitude..
    
    % Update State Vectors (Voltage Angle & Magnitude)
    del(2:nbus) = dTh + del(2:nbus);
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);
            k = k+1;
        end
    end
    Iter = Iter + 1;
    Er = max(abs(M));
end
J;
V
Del = 180/pi*del   % Convert radians to degrees
% Call LoadFlow
%[fb, tb, Pij, Qij] = loadflow(num,V,del,BMva);      
% Call GraphPlot
%graphplot(V,Del,fb,tb,Pij,Qij);

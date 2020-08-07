function [ p_current ] = Get_Intitialp0( G, P )
%GET_INTITIALP0 Summary of this function goes here
%   Detailed explanation goes here

% global L
% global G
L = length(P);

tau_current = G*ones(1,G);

p_next = sdpvar(L, G);
theta_next = sdpvar(L, G);

obj_p = 0;
cons_p = [];

for l = 1:1:L

	for g = 1:1:G
	
		obj_p = obj_p + p_next(l, g);
		
		cons_p = [cons_p, p_next(l, g) >= 0 ]; % Constraint 8d
		
		cons_p = [cons_p, cone( [p_next(l,g), 0.5*(theta_next(l,g)-tau_current(g))], 0.5*(theta_next(l,g)+tau_current(g)) ) ]; % Constraint 40b
		
	end
	
	cons_p = [cons_p, sum(theta_next(l,:)) <= P(l)]; % Constraint 40a
	
end

myops = sdpsettings('solver','sdpt3','verbose',0);

diagnotics = solvesdp(cons_p, -obj_p, myops);
	
OptimalValue_p0 = double(obj_p);

p_current = double(p_next);

theta_current = double(theta_next);

end


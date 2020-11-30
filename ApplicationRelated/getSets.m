function A = getSets(model)


p1 = [1e-4, 1e-6, 1e-8];
p2 = [1e-4, 1e-6, 1e-8];
q2 = [0.2, 0.5, 0.8];

if strcmp(model.name,'log') || strcmp(model.name,'alog')

	A1 = zeros(length(p1),4);
	for i = 1:length(p1)
		fun = @(v)(calculate_actual_value(model, 'math', v, 100, v, 100));
		actual_v = fminsearch(@(v)(abs(log(fun(v)) - log(p1(i)))),10);
		A1(i,:) = [actual_v, 100, actual_v, 100];
	end

	A2 = zeros(length(p2)*length(q2),4);
	counter = 1;
	for i = 1:length(p2)
		for j = 1:length(q2)
			r = Laplace_iCDF(1 - p2(i)/q2(j));
			fun = @(v)(calculate_actual_value(model, 'math', r, 100, -100, v));
			actual_v = fminsearch(@(v)(abs(log(fun(v)) - log(p2(i)))),10);
			A2(counter,:) = [r, 100, -100, actual_v];
			counter = counter + 1;
		end
	end

	A = [A1; A2];

elseif strcmp(model.name,'normal')
    rho = 0.5;
    
	A1 = zeros(length(p1),4);
	for i = 1:length(p1)
		actual_v = MVN_iCDF(p1(i), rho);
		A1(i,:) = [actual_v, 100, actual_v, 100];
	end

	A2 = zeros(length(p2)*length(q2),4);
	counter = 1;
	for i = 1:length(p2)
		for j = 1:length(q2)
			actual_v = MVN_iCDF([p2(i),q2(j)], rho);
			A2(counter,:) = [actual_v(1), 100, -100, actual_v(2)];
			counter = counter + 1;
		end
	end

	A = [A1; A2];

elseif strcmp(model.name,'normal_log_mixture')
    rho = 0.5;
    alpha = 0.5;
    prob = model.tau_cp;
    
	A1 = zeros(length(p1),4);
	for i = 1:length(p1)
		actual_v = MVNLOG_iCDF(p1(i), prob, rho, alpha);
		A1(i,:) = [actual_v, 100, actual_v, 100];
	end

	A2 = zeros(length(p2)*length(q2),4);
	counter = 1;
	for i = 1:length(p2)
		for j = 1:length(q2)
			actual_v = MVNLOG_iCDF([p2(i),q2(j)], prob, rho, alpha);
			A2(counter,:) = [actual_v(1), 100, -100, actual_v(2)];
			counter = counter + 1;
		end
	end

	A = [A1; A2];

end 

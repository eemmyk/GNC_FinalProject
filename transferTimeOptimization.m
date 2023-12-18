function [t_error] = transferTimeOptimization(d_in)
    global theta_f theta_0 tf timeFunction  integralApproximationSteps;
    syms d theta;
    
    %timeFunction = ((6*d*theta - theta^4*((53217828885595963130962570393335*d)/649037107316853453566312041152512 + 3082151914128085151776964925117795/11150372599265311570767859136324180752990208) - theta^2*((102088660786117794935332792382133*d)/20282409603651670423947251286016 + 184767003747656739263195441449317/21778071482940061661655974875633165533184) - theta^6*((3547855259039730875397504692889*d)/1298074214633706907132624082305024 + 205476794275205676785130995007853/22300745198530623141535718272648361505980416) + theta^5*((9515721097982413626635235007929*d)/162259276829213363391578010288128 + 220444152725619773936928989174517/1393796574908163946345982392040522594123776) + d*theta^3 + theta^3*((47578605489912068133176175039645*d)/40564819207303340847894502572032 + 1102220763628098869684644945872585/348449143727040986586495598010130648530944) - theta^4*((34029553595372598311777597460711*d)/81129638414606681695789005144064 + 61589001249218913087731813816439/87112285931760246646623899502532662132736) + 21702051851423/147573952589676412928)/(398600441800000*((- (3547855259039730875397504692889*d)/1298074214633706907132624082305024 - 205476794275205676785130995007853/22300745198530623141535718272648361505980416)*theta^6 + ((9515721097982413626635235007929*d)/162259276829213363391578010288128 + 220444152725619773936928989174517/1393796574908163946345982392040522594123776)*theta^5 + (- (34029553595372598311777597460711*d)/81129638414606681695789005144064 - 61589001249218913087731813816439/87112285931760246646623899502532662132736)*theta^4 + d*theta^3 + 21702051851423/147573952589676412928)^4))^(1/2);
    %timeFunction = ((0.000000000000002508777951886354381858088547658*(theta^3*(1.1729031811226698730080270128789*d + 0.0000000031632184594820754089619961819905) - 1.0*theta^6*(0.0027331682726945404545531099901322*d + 0.000000000009213898120711426515915347087219) - 1.0*theta^4*(0.41944663208613382490466913392257*d + 0.00000000070700706094964492776634202146271) + theta^5*(0.058645159056133493650401350643947*d + 0.00000000015816092297410377044809980909953) + 6.0*d*theta + d*theta^3 - 1.0*theta^2*(5.0333595850336058988560296070709*d + 0.0000000084840847313957391331961042575525) - 1.0*theta^4*(0.081995048180836213636593299703967*d + 0.00000000027641694362134279547746041261657) + 0.00000014705882352941175182300947987812))/(theta^5*(0.058645159056133493650401350643947*d + 0.00000000015816092297410377044809980909953) - 1.0*theta^4*(0.41944663208613382490466913392257*d + 0.00000000070700706094964492776634202146271) - 1.0*theta^6*(0.0027331682726945404545531099901322*d + 0.000000000009213898120711426515915347087219) + d*theta^3 + 0.00000014705882352941175182300947987812)^4)^(1/2);
    %timeFunction = ((6*d*theta - theta^4*((266307113415101805*pi*d)/288230376151711744 - (4906041237903209590593235516275*d)/2535301200456458802993406410752 + 6705668640780052765108028977983285/174224571863520493293247799005065324265472) + d*theta^3 - theta^4*((684465067344555*pi*d)/2251799813685248 - (1234999937906985*d)/5070602400912917605986812821504 + 51704919319757574331031804300025/2722258935367507707706996859454145691648) + theta^5*((13943807851726461*pi*d)/72057594037927936 - (6165113244482127227771269730901*d)/20282409603651670423947251286016 + 421329586785641314854811806081471/43556142965880123323311949751266331066368) - theta^6*((17753807561006787*pi*d)/576460752303423488 - (327069415860213972706215701085*d)/5070602400912917605986812821504 + 447044576052003517673868598532219/348449143727040986586495598010130648530944) + theta^3*((69719039258632305*pi*d)/18014398509481984 - (30825566222410636138856348654505*d)/5070602400912917605986812821504 + 2106647933928206574274059030407355/10889035741470030830827987437816582766592) - theta^2*((2053395202033665*pi*d)/562949953421312 - (3704999813720955*d)/1267650600228229401496703205376 + 155114757959272722993095412900075/680564733841876926926749214863536422912) + 21702051851423/147573952589676412928)/(398600441800000*(((327069415860213972706215701085*d)/5070602400912917605986812821504 - (17753807561006787*d*pi)/576460752303423488 - 447044576052003517673868598532219/348449143727040986586495598010130648530944)*theta^6 + ((13943807851726461*d*pi)/72057594037927936 - (6165113244482127227771269730901*d)/20282409603651670423947251286016 + 421329586785641314854811806081471/43556142965880123323311949751266331066368)*theta^5 + ((1234999937906985*d)/5070602400912917605986812821504 - (684465067344555*d*pi)/2251799813685248 - 51704919319757574331031804300025/2722258935367507707706996859454145691648)*theta^4 + d*theta^3 + 21702051851423/147573952589676412928)^4))^(1/2);
    %timeFunction = (6*d*theta - theta^4*((266307113415101805*pi*d)/288230376151711744 - (4906041237903209590593235516275*d)/2535301200456458802993406410752 + 6705668640780052765108028977983285/174224571863520493293247799005065324265472) + d*theta^3 - theta^4*((684465067344555*pi*d)/2251799813685248 - (1234999937906985*d)/5070602400912917605986812821504 + 51704919319757574331031804300025/2722258935367507707706996859454145691648) + theta^5*((13943807851726461*pi*d)/72057594037927936 - (6165113244482127227771269730901*d)/20282409603651670423947251286016 + 421329586785641314854811806081471/43556142965880123323311949751266331066368) - theta^6*((17753807561006787*pi*d)/576460752303423488 - (327069415860213972706215701085*d)/5070602400912917605986812821504 + 447044576052003517673868598532219/348449143727040986586495598010130648530944) + theta^3*((69719039258632305*pi*d)/18014398509481984 - (30825566222410636138856348654505*d)/5070602400912917605986812821504 + 2106647933928206574274059030407355/10889035741470030830827987437816582766592) - theta^2*((2053395202033665*pi*d)/562949953421312 - (3704999813720955*d)/1267650600228229401496703205376 + 155114757959272722993095412900075/680564733841876926926749214863536422912) + 21702051851423/147573952589676412928)/(398600441800000*(((327069415860213972706215701085*d)/5070602400912917605986812821504 - (17753807561006787*d*pi)/576460752303423488 - 447044576052003517673868598532219/348449143727040986586495598010130648530944)*theta^6 + ((13943807851726461*d*pi)/72057594037927936 - (6165113244482127227771269730901*d)/20282409603651670423947251286016 + 421329586785641314854811806081471/43556142965880123323311949751266331066368)*theta^5 + ((1234999937906985*d)/5070602400912917605986812821504 - (684465067344555*d*pi)/2251799813685248 - 51704919319757574331031804300025/2722258935367507707706996859454145691648)*theta^4 + d*theta^3 + 21702051851423/147573952589676412928)^4);

    timeFunction_n = subs(timeFunction, d, d_in);
    %transferTime = @(angle) double(subs(timeFunction_n, theta, angle));


    %Transfer Time
    theta_vec = linspace(theta_0, theta_f, integralApproximationSteps);
    timeCurve = zeros(2, integralApproximationSteps);
    for i = 1:integralApproximationSteps
        timeCurve(:,i) = [theta_vec(i); double(subs(timeFunction_n, theta, theta_vec(i)))];
    end
    
    time_t = abs(trapz(timeCurve(1, :), timeCurve(2,:)));

    %time_t = integral(transferTime, theta_0, theta_f, "AbsTol",1e-2, "RelTol", 1e-2);
    t_error = time_t - tf;

    fprintf("guessed d: %.3f, remaining error: %.2f\n", d_in, t_error);

end


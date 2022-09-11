#define _CRT_SECURE_NO_WARNINGS

#include "integrator.h"
#include "integrator_dp54.h"
#include <functional>

IMathModel* test_runge() {
	//gauss G(-1, 4);
	vector<double> m = { 10, 10, 10, 3 };

	//MotionModel* model = new MotionModel(m, &G);

	attractor* modelA = new attractor(m);
	/*eylerArgs args1(modelA, 10, 0.01, 0.01, 0);
	eyler c1(args1); */

	/*eylerCorrectionArgs args2(model, 100, 0.01, 1, 0);
	eylerCorrection c1(args2); */

	rungeArgs args(modelA, 1000, 0.01);
	runge c1(args);

	c1.integration();
	//cin.get();

	return modelA;
}

IMathModel* test_dp54() {
	double 
		t0 = 0, 
		tf = 100, 
		h0 = 0.001,
		epsMax = 0.0000001,
		k = 1;
	vector<double> m = { 0.994, 0, 0, -2.00158510637908252240537862224  };
	IMathModel* modelA = new dp54_Model(m);
	dp54_args args(modelA, t0, tf, 0, 0, h0, epsMax, k);
	dp54_integrator c1(args);
	c1.integration();	
	return modelA;
}

IMathModel* test_dp54_pendulum() {
	double
		t0 = 0,
		tf = 160,
		h0 = 0.001,
		epsMax = 0.00001,
		k = 0.01;
	//IMathModel* modelA = new pendulum_ideal_Model(t0);    //1
	//IMathModel* modelA = new pendulum_ideal_spring_Model(t0);    //2
	//IMathModel* modelA = new pendulum_real_Model(t0);     //3
	//IMathModel* modelA = new pendulum_real_spring_slide_Model(t0);    //4
	IMathModel* modelA = new pendulum_real_viscous_Model(t0);   //5
	dp54_args args(modelA, t0, tf, 0, 0, h0, epsMax, k);
	dp54_integrator c1(args); // 1 2  for ideal
	//dp54_integrator_2 c1(args); // 3 4 5 for real
	c1.integration();
	return modelA;
}

void attractor_show(IMathModel* model, function<vector<double> (vector<double>)> map) {
	attractor_show_args args_attr;
	args_attr.example_to_show = false;
	args_attr.example_continiously_draw = true;
	auto v = model->X_all();
	cout << "size : " << v.size() << endl;
	vector<vector<double>> points(v.size());
	std::transform(v.begin(), v.end(), points.begin(), map);
	args_attr.points = points;
	args_attr.print_dots = false;
	attractor_show(args_attr);
}


void main() {
	srand(time(nullptr));
	//attractor_show(test_runge()); // todo coefficiens * 40 to remove
	
	// for DP
	/*attractor_show(test_dp54(), [](vector<double> v) -> vector<double> {
		return vector<double>{ v[0] * 40, 0, v[1] * 40 };
	});*/


	// for pendulum
	
	attractor_show(test_dp54_pendulum(), [](vector<double> v) -> vector<double> {
		//return vector<double>{ v[0] * 3, 0, v[1] * 3 };                             // 1 2 3 4 
		return vector<double>{ v[0]*0.8, 0, v[1]*0.1 };                         // for pendulum_real_viscous_Model
});
	std::cin.get();
}
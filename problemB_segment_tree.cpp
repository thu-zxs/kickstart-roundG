/* 
 * authored by github.com/thu-zxs after the competition.
 */
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
using namespace std;
template <typename T>
class Interval
{
	public:
		pair<T, T> range;
		pair<bool, bool> closure; 
		Interval(): range(pair<T, T>(-1, -1)), closure(pair<bool,bool>(0,0)) {}
		Interval(pair<T, T> &range_tmp, pair<bool,bool>& closure_tmp):
			range(range_tmp), closure(closure_tmp) {}
		Interval(T f, T s, bool cf, bool cs) {
			range = make_pair(f, s);
			closure = make_pair(cf, cs);
		}
		bool intersect(Interval& sec_intv) {
			if ( !(
					(range.first > sec_intv.range.second || range.second < sec_intv.range.first) 
					|| 
					(
					 	(range.first == sec_intv.range.second&&!(closure.first==true&&sec_intv.closure.second==true))
						 ||
						(range.second==sec_intv.range.first&&!(closure.second==true&&sec_intv.closure.first==true))
					)
				  ))
				return true;
			else
				return false;
		}
		bool contains(Interval& sec_intv) {
			if (
					(range.first < sec_intv.range.first||range.first==sec_intv.range.first&&closure.first>=sec_intv.closure.first) && (range.second > sec_intv.range.second||range.second==sec_intv.range.second&&closure.second>=sec_intv.closure.second)
				)
				return true;
			else
				return false;
		}
		void print() {
			if(closure.first) cout<<'[';
			else cout<<'(';
			cout<<range.first<<','<<range.second;
			if(closure.second) cout<<']';
			else cout<<')';
		}
};

template <typename T>
class TNode
{
	public:
		Interval<T> interval;
		vector<Interval<T> > subset;
		TNode(Interval<T> &intv) {
			interval.range = intv.range;
			interval.closure = intv.closure;
		}
		TNode() {} // trivial node
};

template <typename T1> class SegmentTree
{ 
	public:
		vector<TNode<T1> > nodes; // store tree structure
		vector<Interval<T1> > ele_int; // store leaf node.
		vector<T1> cnt;

		SegmentTree(vector<T1> &eles) {
			Init(eles);
		}

		void printTree() {
			for(auto &n: nodes) {
				n.interval.print();
				cout << ": ";
				for (auto &intv: n.subset) {
					intv.print();
					cout << ' ';
				}
				cout << endl;
			}
		}

		void printCnt() {
			for(int i = 0; i < ele_int.size(); i++) {
				ele_int[i].print();
				cout << ": ";
				cout << cnt[i] << endl;
			}
		}

		void insert(T1 root, Interval<T1> &intv) {
			// if( nodes[root].start == -1) return;
			// intv.print(); cout << endl;
			// cout << root << ' '; nodes[root].interval.print(); cout << ' '; nodes[2*root].interval.print(); cout << ' '; nodes[2*root+1].interval.print(); cout << ' '; intv.print(); cout << endl;
			if( intv.contains(nodes[root].interval) ) {
				nodes[root].subset.push_back(intv);
				// cout << root << ' '; nodes[root].interval.print(); cout << ' '; nodes[2*root].interval.print(); cout << ' '; nodes[2*root+1].interval.print(); cout << ' '; intv.print(); cout << endl;
			} else {
				// cout << root << ' '; nodes[root].interval.print(); cout << ' '; nodes[2*root].interval.print(); cout << ' '; nodes[2*root+1].interval.print(); cout << ' '; intv.print(); cout << endl;
				if (nodes[2*root].interval.intersect(intv)) insert(2*root, intv);
				if (nodes[2*root+1].interval.intersect(intv)) insert(2*root+1, intv);
			}
		}

		void query_ele_int() {
			for(auto &ele: ele_int) {
				if(nodes[1].interval.contains(ele)) cnt.push_back(query_each_ele_int(1, ele));
				else cnt.push_back(0);
			}
		}

		T1 query_each_ele_int(T1 root, Interval<T1> &ele) {
			// if (nodes[root].start == -1) return 0;
			T1 cnt = (T1)nodes[root].subset.size();
			if ( 2*root < nodes.size() && nodes[2*root].interval.contains(ele) )
				cnt += query_each_ele_int(2*root, ele);
			if ( 2*root+1 < nodes.size() && nodes[2*root+1].interval.contains(ele) )
				cnt += query_each_ele_int(2*root+1, ele);
			return cnt;
		}

	private:
		void Init(vector<T1> &eles) { // `eles` is sorted
			int n = eles.size();
			for(int i = 0; i < n; i++) {
				Interval<T1> tmp(eles[i], eles[i], true, true); 
				ele_int.push_back( tmp );
				if(i<n-1) {
					Interval<T1> tmp1(eles[i], eles[i+1], false, false);
					ele_int.push_back( tmp1 );
				}
			}
			int size_ele_int = ele_int.size();
			nodes.resize(size_ele_int * 4); // for array storage of binary tree.
			cout << "nodes size: " << size_ele_int * 4 << endl;
			build(ele_int, 1, 0, size_ele_int-1);
		}

		void build(vector<Interval<T1> > &ele_int, T1 root, int ind_start, int ind_end) { // by default settings, root is at position 1. O(N)
			if(ind_start == ind_end) {
				nodes[root].interval.range.first = ele_int[ind_start].range.first;
				nodes[root].interval.range.second = ele_int[ind_end].range.second;
				nodes[root].interval.closure.first = ele_int[ind_start].closure.first;
				nodes[root].interval.closure.second = ele_int[ind_end].closure.second;
			} else {
				int ind_mid = (ind_start+ind_end) / 2;
				build(ele_int, 2*root, ind_start, ind_mid);
				build(ele_int, 2*root+1, ind_mid+1, ind_end);
				nodes[root].interval.range.first = nodes[2*root].interval.range.first;
				nodes[root].interval.range.second = nodes[2*root+1].interval.range.second;
				nodes[root].interval.closure.first = nodes[2*root].interval.closure.first;
				nodes[root].interval.closure.second = nodes[2*root+1].interval.closure.second;
			}
		}
};

template <typename T>
int binary_search(vector<T> &array, T k, int left, int right) { // array in descend order
	// cout << array[0] << endl;
	if (array[0] < k) return -1;

	int mid = (left+right) / 2;

	// cout << mid << ' ' << array[mid] << endl;

	if ( mid < array.size()-1 )
		{if (array[mid] >= k && array[mid+1] < k) return mid;}
	else if (mid == array.size() - 1)
		return mid;

	if (left >= right) return -1;

	if (array[mid] < k)
		return binary_search(array, k, left, mid);
	else if (array[mid] >= k)
		return binary_search(array, k, mid+1, right);
}

template <typename T>
T count(int pos, vector<T> &L, vector<T> &R) {
	T cnt = 0;
	for (int i = 0; i < L.size(); i++) {
		if (pos < L[i]) cnt += R[i]-L[i]+1;
		else if(pos >= L[i] && pos <= R[i]) cnt += R[i]-pos+1;
	}
	return cnt;
}

template <typename T>
T binary_search(T k, T left, T right, vector<T> &L, vector<T> &R) { // array in descend order

	T mid = (left+right) / 2;
	T cnt = count(mid, L, R);

	// cout << mid << ' ' << array[mid] << endl;

	// if ( mid < array.size()-1 )
	// 	{if (array[mid] >= k && array[mid+1] < k) return mid;}
	// else if (mid == array.size() - 1)
	// 	return mid;
	//
	if ( cnt >= k && count(mid+1, L, R) < k) return mid;

	if (left >= right) return -1;

	if (cnt < k)
		return binary_search(k, left, mid, L, R);
	else if (cnt >= k)
		return binary_search(k, mid+1, right, L, R);
}

/*
int main() {
	ifstream fin("inb.txt");
	ofstream fout("outb.txt");
	int T;
	fin >> T;
	for (int t = 0; t < T; t++) {
		int N;
		int Q;
		long long x1,x2,a1,b1,c1,m1,  y1,y2,a2,b2,c2,m2,  z1,z2,a3,b3,c3,m3;
		fin >> N >> Q;
		fin >>x1>>x2>>a1>>b1>>c1>>m1;  
		fin >>y1>>y2>>a2>>b2>>c2>>m2; 
		fin >>z1>>z2>>a3>>b3>>c3>>m3;

		vector<long long> X(N, 0), Y(N, 0), Z(Q, 0), L(N, 0), R(N, 0);
		vector<long long> K(Q, 0);
		X[0] = x1; X[1] = x2;
		Y[0] = y1; Y[1] = y2;
		Z[0] = z1; Z[1] = z2;
		for (int i=2; i < N; i++) {
			X[i] = (a1 * X[i-1] + b1 * X[i-2] + c1) % m1;
			Y[i] = (a2 * Y[i-1] + b2 * Y[i-2] + c2) % m2;
		}
		for (int i=2; i < Q; i++) {
			Z[i] = (a3 * Z[i-1] + b3 * Z[i-2] + c3) % m3;
		}

		for (int i=0; i < N; i++) {
			L[i] = min(X[i], Y[i]) + 1;
			R[i] = max(X[i], Y[i]) + 1;
		}
		for (int i=0; i < Q; i++) {
			K[i] = Z[i] + 1;
		}
		auto minLp = min_element(L.begin(), L.end());
		auto maxRp = max_element(R.begin(), R.end());
		auto minL = *minLp;
		auto maxR = *maxRp;

		long long ans = 0;
		long long all_count = count(minL, L, R);
		for (int k = 0; k < Q; k++) {
			long long tmp;
			if (all_count < K[k]) tmp = -1;
			else tmp = binary_search(K[k], minL, maxR, L, R);
			if(tmp == -1) continue;
			ans += (k+1)*tmp;
		}
		fout << "Case #" << t+1 << ": " << ans << endl;
	}
	return 0;
}
*/
int main() {

	ifstream fin("inb.txt");
	ofstream fout("outb_seg.txt");
	int T;
	fin >> T;
	for (int t = 0; t < T; t++) {
		int N;
		int Q;
		long long x1,x2,a1,b1,c1,m1,  y1,y2,a2,b2,c2,m2,  z1,z2,a3,b3,c3,m3;
		fin >> N >> Q;
		fin >>x1>>x2>>a1>>b1>>c1>>m1;  
		fin >>y1>>y2>>a2>>b2>>c2>>m2; 
		fin >>z1>>z2>>a3>>b3>>c3>>m3;

		vector<long long> X(N, 0), Y(N, 0), Z(Q, 0), L(N, 0), R(N, 0);
		vector<long long> K(Q, 0);
		X[0] = x1; X[1] = x2;
		Y[0] = y1; Y[1] = y2;
		Z[0] = z1; Z[1] = z2;
		for (int i=2; i < N; i++) {
			X[i] = (a1 * X[i-1] + b1 * X[i-2] + c1) % m1;
			Y[i] = (a2 * Y[i-1] + b2 * Y[i-2] + c2) % m2;
		}
		for (int i=2; i < Q; i++) {
			Z[i] = (a3 * Z[i-1] + b3 * Z[i-2] + c3) % m3;
		}

		for (int i=0; i < N; i++) {
			L[i] = min(X[i], Y[i]) + 1;
			R[i] = max(X[i], Y[i]) + 1;
		}
		for (int i=0; i < Q; i++) {
			K[i] = Z[i] + 1;
		}

		vector<Interval<long long> > intervals;
		intervals.resize(L.size());
		for (int i = 0; i < L.size(); i++) {
			// cout << L[i] << endl;
			intervals[i].range.first = L[i];
			intervals[i].range.second = R[i];
			intervals[i].closure.first = true;
			intervals[i].closure.second = true;
		}

		set<long long> tmp;
		vector<long long> eles;
		for (int i = 0; i < L.size(); i++) {
			tmp.insert(L[i]);
			tmp.insert(R[i]);
		}
		eles.assign(tmp.begin(), tmp.end());

		SegmentTree<long long> tree(eles);
		for (int i = 0; i < intervals.size(); i++) {
			tree.insert(1, intervals[i]);
		}
		// cout << t << endl;
		// tree.printTree();
		tree.query_ele_int();
		// tree.printCnt();
		int n = tree.cnt.size();
		vector<long long> suffix_sum(n, 0);
		suffix_sum[n-1] = tree.cnt[n-1];
		for(int i = n-2; i >= 0; i--) {
			int range_len = tree.ele_int[i].range.second - tree.ele_int[i].range.first;
			if (range_len == 0)
				suffix_sum[i] = tree.cnt[i] + suffix_sum[i+1];
			else if (range_len > 1)
				suffix_sum[i] = tree.cnt[i]*(range_len-1) + suffix_sum[i+1];
			else
				suffix_sum[i] = suffix_sum[i+1];
			// cout << suffix_sum[i] << endl;
		}

		long long ans = 0;
		for (int i = 0; i < K.size(); i++) {
			int pre_ans = binary_search(suffix_sum, K[i], 0, suffix_sum.size()-1);
			if (pre_ans==-1) {
				ans += 0;
			}
			else {
				// tree.ele_int[pre_ans].print(); cout << ' ' << K[i] << ' ' << suffix_sum[pre_ans] << endl;
				long long range_len = tree.ele_int[pre_ans].range.second - tree.ele_int[pre_ans].range.first;
				if (range_len == 0) ans += (i+1)*tree.ele_int[pre_ans].range.first;
				else if (range_len > 1) {
					ans += (i+1)*(tree.ele_int[pre_ans].range.first+(suffix_sum[pre_ans]-K[i])/tree.cnt[pre_ans] + 1);
				}
			}
		}
		fout << "Case #" << t+1 << ": " << ans << endl;
	}
	return 0;
}

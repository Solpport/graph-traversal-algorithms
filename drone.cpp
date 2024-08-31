// Project Identifier: 1761414855B69983BD8035097EFBD312EB0527F0
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <queue>

struct point
{
	double x;
	double y;
};

struct points
{
	std::vector<point> pts;
	std::vector<std::size_t> medical_campus; // QIII
	std::vector<std::size_t> main_campus; // QI, QII, QIV
	std::vector<std::size_t> border; // origin + negative x/y axes
};

points get_points();
void do_mst(const points& pts);
void do_fasttsp(const points& pts);
void do_opttsp(const points& pts);

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "Invalid command line.\n";
		return 1;
	}

	enum class mode_t { MST, FASTTSP, OPTTSP } mode;

	std::string_view arg1 = argv[1];
	if (arg1 == "-m" || arg1 == "--mode")
	{
		if (argc != 3)
		{
			std::cerr << "Invalid command line.\n";
			return 1;
		}

		std::string_view arg2 = argv[2];
		if (arg2 == "MST")
			mode = mode_t::MST;
		else if (arg2 == "FASTTSP")
			mode = mode_t::FASTTSP;
		else if (arg2 == "OPTTSP")
			mode = mode_t::OPTTSP;
		else
		{
			std::cerr << "Invalid command line.\n";
			return 1;
		}
	}
	else if (arg1 == "-h" || arg1 == "--help")
	{
		std::cout << "Help menu\n";
		return 0;
	}
	else
	{
		std::cerr << "Invalid command line.\n";
		return 1;
	}

	std::cout << std::fixed << std::setprecision(2);

	try
	{
		points pts = get_points();
		switch (mode)
		{
		case mode_t::MST:
			do_mst(pts);
			return 0;
		case mode_t::FASTTSP:
			do_fasttsp(pts);
			return 0;
		case mode_t::OPTTSP:
			do_opttsp(pts);
			return 0;
		}
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
	
}

points get_points()
{
	std::size_t len;
	if (!(std::cin >> len))
		throw std::runtime_error("Invalid input");

	points pts;
	pts.pts.resize(len);
	pts.medical_campus.reserve(len);
	pts.main_campus.reserve(len);
	pts.border.reserve(len);
	for (std::size_t i = 0; i < len; ++i)
	{
		point p;
		if (!(std::cin >> p.x >> p.y))
			throw std::runtime_error("Invalid input");
		
		pts.pts[i] = p;
		
		if (p.y > 0) // QI, QII
		{
			pts.main_campus.push_back(i);
		}
		else if (p.y < 0)
		{
			if (p.x > 0) // QIV
				pts.main_campus.push_back(i);
			else if (p.x < 0) // QIII
				pts.medical_campus.push_back(i);
			else // negative y axis
				pts.border.push_back(i);
		}
		else
		{
			if (p.x <= 0) // negative x axis or origin
				pts.border.push_back(i);
			else
				pts.main_campus.push_back(i);
		}
	}
	
	return pts;
}

double dist(const points &pts, unsigned int from, unsigned int to)
{
	double xdiff = pts.pts[to].x - pts.pts[from].x;
	double ydiff = pts.pts[to].y - pts.pts[from].y;
	return std::sqrt(xdiff * xdiff + ydiff * ydiff);
}

struct prims_graph
{
	struct edge
	{
		unsigned int from, to;
	};

	prims_graph(const points &p) : graph(p.pts.size()), pts(p)
	{
		// edges.reserve((((pts.main_campus.size() - 1) * pts.main_campus.size() + 1) / 2) + (pts.main_campus.size() * pts.border.size()) +
		// 			  (((pts.medical_campus.size() - 1) * pts.medical_campus.size() + 1) / 2) + (pts.medical_campus.size() * pts.border.size()));
		auto setup_zone = [&](const std::vector<std::size_t> &zone)
		{
			if (zone.empty())
				return;

			for (std::size_t i = 0; i < zone.size(); ++i)
			{
				auto &cur = graph[zone[i]];
				std::size_t j = 0;
				for (; j < i; ++j)
					cur.push_back((unsigned int)zone[j]);
				for (++j; j < zone.size(); ++j)
					cur.push_back((unsigned int)zone[j]);

				for (j = 0; j < pts.border.size(); ++j)
				{
					cur.push_back((unsigned int)pts.border[j]);
					graph[pts.border[j]].push_back((unsigned int)zone[i]);
				}
			}
		};

		setup_zone(pts.main_campus);
		setup_zone(pts.medical_campus);
	}

	// [path, distance]
	std::pair<std::vector<edge>, double> prims()
	{
		std::pair<std::vector<edge>, double> res;
		res.first.reserve(pts.pts.size());
		res.second = 0;

		std::vector<bool> included(pts.pts.size());
		included[0] = true;

		auto cmp = [this](const edge &a, const edge &b) { return dist(pts, a.from, a.to) > dist(pts, b.from, b.to); };
		std::vector<edge> buff(graph.front().size());
		for (std::size_t i = 0; i < buff.size(); ++i)
		{
			buff[i].from = 0;
			buff[i].to = graph.front()[i];
		}
		std::priority_queue<edge, std::vector<edge>, decltype(cmp)> queue(cmp, std::move(buff));
		while (!queue.empty())
		{
			edge e = queue.top();
			queue.pop();
			if (!included[e.to])
			{
				included[e.to] = true;
				res.first.push_back(e);
				res.second += dist(pts, e.from, e.to);
				for (auto next : graph[e.to])
					if (!included[next])
						queue.push({e.to, next});
			}
		}

		if (res.first.size() != pts.pts.size() - 1)
			res.first.clear();

		return res;
	}

private:
	std::vector<std::vector<unsigned int>> graph;
	const points &pts;
};

void do_mst(const points& pts)
{
	prims_graph graph(pts);
	auto res = graph.prims();

	if (res.first.empty())
	{
		std::cerr << "Cannot construct MST\n";
		std::exit(1);
	}

	std::cout << res.second << '\n';
	for (auto e : res.first)
		if (e.from < e.to)
			std::cout << e.from << ' ' << e.to << '\n';
		else
			std::cout << e.to << ' ' << e.from << '\n';
}

std::pair<std::vector<unsigned int>, double> nearest_neigbhor(const points &pts)
{
	std::pair<std::vector<unsigned int> , double> res{};
	if (pts.pts.empty())
		return res;
	
	std::vector<unsigned int> unvisited(pts.pts.size() - 1);
	for (unsigned int i = 1; i < pts.pts.size(); ++i)
		unvisited[i - 1] = i;
	
	res.first.reserve(pts.pts.size());
	res.first.push_back(0);

	while (unvisited.size() != 0)
	{
		double min_dist = std::numeric_limits<double>::infinity();
		unsigned int min_i = 0;

		unsigned int prev = res.first.back();
		for (std::size_t i = 0; i < unvisited.size(); ++i)
		{	
			double d = dist(pts, prev, unvisited[i]);
			if (d < min_dist)
			{
				min_dist = d;
				min_i = static_cast<unsigned int>(i);
			}
		}

		res.first.push_back(unvisited[min_i]);
		res.second += min_dist;
		std::swap(unvisited[min_i], unvisited[unvisited.size() - 1]);
		unvisited.pop_back();
	}

	res.second += dist(pts, res.first.back(), 0);

	return res;
}

double delta_segment(const std::vector<unsigned int> &path, std::size_t i, std::size_t j, const points &pts)
{
	double d_old = dist(pts, path[i - 1], path[i]) + dist(pts, path[i], path[i + 1]) +
				   dist(pts, path[j - 1], path[j]) + dist(pts, path[j], path[(j + 1) % path.size()]);

	double d_new = dist(pts, path[i - 1], path[j]) + dist(pts, path[j], path[j - 1]) +
				   dist(pts, path[i + 1], path[i]) + dist(pts, path[i], path[(j + 1) % path.size()]);

	return d_new - d_old;
}

void two_opt(std::vector<unsigned int> &path, double cur_dist, const points &pts)
{
	if (path.size() <= 2)
		return;
	
	bool improved;
	do
	{
		improved = false;
		for (std::size_t i = 1; i < path.size() - 2; ++i)
			for (std::size_t j = i + 1; j < path.size(); ++j)
			{
				double delta = delta_segment(path, i, j, pts); // calculate delta
				if (delta < 0)
				{
					std::reverse(path.begin() + i, path.begin() + j + 1);; // reverse segnment
					cur_dist += delta;
					improved = true;
				}
			}
	} while (improved);
}

void do_fasttsp(const points& pts)
{
	auto initial = nearest_neigbhor(pts);
	
	two_opt(initial.first, initial.second, pts);
	
	std::cout << initial.second << '\n';
	for (auto i : initial.first)
		std::cout << i << ' ';
	std::cout << '\n';
}

struct perm_generator
{
	perm_generator(const points& pts) : best_dist(std::numeric_limits<double>::infinity()), pts(pts)
	{	
		std::vector<unsigned int> path(pts.pts.size());
		for (std::size_t i = 0; i < path.size(); ++i)
			path[i] = static_cast<unsigned int>(i);
		gen(path, 1, 0);
	}

	std::vector<unsigned int> best;
	double best_dist;

private:
	double mst_cost(const std::vector<unsigned int> &path, std::size_t len) const
	{
		if (len >= path.size())
			return 0.0; // No unvisited nodes left.

		double total_mst_cost = 0.0;
		std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>, std::greater<>> min_heap;
		std::vector<bool> visited(path.size(), false);
		std::vector<double> min_edge_to_tree(path.size(), std::numeric_limits<double>::infinity());

		// Start MST with the last node in the current path
		unsigned int start_node = path[len - 1];
		min_edge_to_tree[start_node] = 0;
		min_heap.push({0.0, start_node});

		while (!min_heap.empty())
		{
			auto [cost, u] = min_heap.top();
			min_heap.pop();

			if (visited[u])
				continue;
			visited[u] = true;
			total_mst_cost += cost;

			// Update the minimum edge cost to each unvisited node
			for (std::size_t i = len; i < path.size(); ++i)
			{
				unsigned int v = path[i];
				if (!visited[v] && dist(pts, u, v) < min_edge_to_tree[v])
				{
					min_edge_to_tree[v] = dist(pts, u, v);
					min_heap.push({min_edge_to_tree[v], v});
				}
			}
		}

		return total_mst_cost;
	}

	bool promising(const std::vector<unsigned int> &path, std::size_t len, double cur_dist) const
	{
		double min = mst_cost(path, len);
		return cur_dist + min < best_dist;
	}

	void gen(std::vector<unsigned int>& path, std::size_t len, double cur_dist)
	{
		if (len == path.size())
		{
			cur_dist += dist(pts, path[len - 1], path[0]); // return to start
			if (cur_dist < best_dist)
			{
				best_dist = cur_dist;
				best = path;
			}

			return;
		}

		if (!promising(path, len, cur_dist))
			return;

		for (std::size_t i = len; i < path.size(); ++i)
		{
			std::swap(path[len], path[i]);
			gen(path, len + 1, cur_dist + dist(pts, path[len - 1], path[len]));
			std::swap(path[len], path[i]);
		}
	}
	
	const points &pts;
};

void do_opttsp(const points& pts)
{
	perm_generator gen(pts);
	std::cout << gen.best_dist << '\n';
	for (auto i : gen.best)
		std::cout << i << ' ';
	std::cout << '\n';
}

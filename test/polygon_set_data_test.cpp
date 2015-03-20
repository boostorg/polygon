// Boost.Polygon library polygon_set_data_test.cpp file

//          Copyright Andrii Sydorchuk 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#define BOOST_TEST_MODULE POLYGON_SET_DATA_TEST
#include <vector>

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>

#include "boost/polygon/polygon.hpp"
using namespace boost::polygon;

typedef boost::mpl::list<int> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(polygon_set_data_test, T, test_types) {
    typedef point_data<T> point_type;
    typedef polygon_data<T> polygon_type;
    typedef polygon_with_holes_data<T> polygon_with_holes_type;
    typedef polygon_set_data<T> polygon_set_type;

    polygon_set_type pset;
    std::vector<point_type> outbox;
    outbox.push_back(point_type(0, 0));
    outbox.push_back(point_type(100, 0));
    outbox.push_back(point_type(100, 100));
    outbox.push_back(point_type(0, 100));
    pset.insert_vertex_sequence(outbox.begin(), outbox.end(), COUNTERCLOCKWISE, false);
    std::vector<point_type> inbox;
    inbox.push_back(point_type(20, 20));
    inbox.push_back(point_type(80, 20));
    inbox.push_back(point_type(80, 80));
    inbox.push_back(point_type(20, 80));
    pset.insert_vertex_sequence(inbox.begin(), inbox.end(), COUNTERCLOCKWISE, true);

    BOOST_CHECK(!pset.empty());
    BOOST_CHECK(!pset.sorted());
    BOOST_CHECK(pset.dirty());
    BOOST_CHECK_EQUAL(8, pset.size());

    std::vector<polygon_with_holes_type> vpoly;
    pset.get(vpoly);
    BOOST_CHECK_EQUAL(1, vpoly.size());

    polygon_with_holes_type poly = vpoly[0];
    BOOST_CHECK_EQUAL(5, poly.size());
    BOOST_CHECK_EQUAL(1, poly.size_holes());
}

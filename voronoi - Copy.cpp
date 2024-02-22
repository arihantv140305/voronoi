#include <vector>
#include <algorithm>
#include <cmath>
#include<iostream>
using namespace std;

#include<bits/stdc++.h>

double directrix = 0.0; 
double epsilon = 1e-10; // or whatever small value you consider "close enough"
struct Event;
struct Point {
    double x;
    double y;
    Event* circle_event;
    Point(){};
    Point(double x, double y) : x(x), y(y) {}
    bool operator==(const Point& other) const {
        return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
    }
};


struct breakpoint {
    public:
    Point f1, f2;

    breakpoint(Point f1, Point f2) : f1(f1), f2(f2) {
        calculate_intersection();
    }
    breakpoint(){};
    double calculate_intersection() const {
        double x1 = f1.x, y1 = f1.y - directrix, x2 = f2.x, y2 = f2.y - directrix;

        double a = y2 - y1;
        if(a==0) {
            return (x1+x2)/2;
        }
        double b = -2*(x1*y2 - x2*y1);
        double c = y2*x1*x1 - y1*x2*x2 - y1*y2*(y2 - y1);     
        double d = b*b - 4*a*c;
        return (-b + sqrt(d))/(2*a);
    }
    bool operator<(const breakpoint& other) const {
        return calculate_intersection() < other.calculate_intersection();
    }
    bool operator==(const breakpoint& other) const {
        return calculate_intersection() == other.calculate_intersection();
    }
    void printit(){
        cout<<"("<<f1.x<<" "<<f1.y<<")"<<" ("<<f2.x<<" "<<f2.y<<")"<<endl;
    }
};

Point getCircleEvent(breakpoint l, breakpoint r){
    if(!(l.f2 == r.f1)){
        l.printit();
        r.printit();
        exit(0);
    }
    Point p1 = l.f1, p2 = l.f2, p3 = r.f2;
    // Find lowest point of circumcircle of p1, p2, p3
    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;
    double x3 = p3.x, y3 = p3.y;

    double A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
    double B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
    double C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);

    double centerX = -B / (2 * A);
    double centerY = -C / (2 * A);

    double radius = sqrt((x1 - centerX) * (x1 - centerX) + (y1 - centerY) * (y1 - centerY));
    return Point(centerX, centerY - radius);
}

struct Event {
    Point p;
    breakpoint b1;
    breakpoint b2;
    bool is_site_event;
    Event(Point p, bool is_site_event) : p(p), is_site_event(is_site_event) {}
    bool operator<(const Event& other) const {
        return p.y > other.p.y;
    }
    bool operator==(const Event& other) const {
        return p == other.p;
    }
    Event(breakpoint b1, breakpoint b2, bool is_site_event) : b1(b1), b2(b2), is_site_event(is_site_event) {
        p = getCircleEvent(b1, b2);
    }
};

bool check_clockwise(breakpoint bk1, breakpoint bk2){
    Point a = bk1.f1;
    Point b = bk1.f2;
    Point c = bk2.f2;
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x) < 0;
}

class Beachline{
    public:
    // struct CompareBreakpoints{
    //     bool operator()(const breakpoint& a,const breakpoint& b) {
    //         return a.calculate_intersection() < b.calculate_intersection();
    //     }
    // };
    multiset<breakpoint> beachline;
    multiset<Event> events;
    void handle_site_event(Point p0) {
        if (beachline.empty()) {
            // If the beachline is empty, just insert p0 as a breakpoint with itself
            breakpoint k = breakpoint(p0, p0);  
            beachline.insert(k);
        } else if (beachline.size() == 1) {
            // If there is only one point in the beachline, insert two new breakpoints
            auto bkpt = beachline.begin();
            beachline.insert(breakpoint(bkpt->f1, p0));
            beachline.insert(breakpoint(p0, bkpt->f2));
            beachline.erase(bkpt);
        } else {
            // If there are at least two points in the beachline
            breakpoint k = breakpoint(p0, p0);
            auto it = beachline.lower_bound(k); // Find the first breakpoint that does not compare less than p0 
            auto bkpt1 = std::prev(it);
            auto bkpt2 = it;
            if (it == beachline.begin()){ 
                bkpt1 = bkpt2;
                breakpoint nbp1 = breakpoint(p0, bkpt1->f1);
                if(getCircleEvent(nbp1, *bkpt1).y < directrix && check_clockwise(nbp1, *bkpt1)){
                    cout<<"Circle inserted at leftmost end, using\n";
                    auto x = *bkpt1;
                    x.printit();
                    nbp1.printit();
                    Event new_event = Event(nbp1, *bkpt1, false);
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }
                beachline.insert(breakpoint(bkpt1->f1, p0));
                beachline.insert(breakpoint(p0, bkpt1->f1));
            }
            else if(it == beachline.end()) {
                bkpt2 = bkpt1;
                breakpoint nbp2 = breakpoint(bkpt2->f2, p0);
                if(getCircleEvent(*bkpt2, nbp2).y < directrix && check_clockwise(*bkpt2, nbp2)){
                    cout<<"Circle inserted at rightmost\n";
                    Event new_event = Event(*bkpt2, nbp2, false);
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }

                beachline.insert(breakpoint(bkpt1->f2, p0));
                beachline.insert(breakpoint(p0, bkpt1->f2));
            }
            else{
                if(bkpt1->f2 == bkpt2->f1 && bkpt1->f1 == bkpt2->f2){
                    cout<<"NO CIRCLE ABOVE"<<endl;
                }
                else{
                    cout<<"CIRCLE MAY BE ABOVE"<<endl;   
                    Event e = Event(*bkpt1, *bkpt2, false);
                    if(events.find(e) != events.end()){
                        cout<<"Circle found and erased ("<<e.p.x<<" "<<e.p.y<<")"<<endl;
                        events.erase(e);
                    }
                }
                breakpoint nbp1 = breakpoint(bkpt1->f2, p0);
                breakpoint nbp2 = breakpoint(p0, bkpt2->f1);
                if(getCircleEvent(*bkpt1, nbp1).y < directrix && check_clockwise(*bkpt1, nbp1)){
                    cout<<"Circle inserted using\n";
                    auto x = *bkpt1;
                    x.printit();
                    nbp1.printit();
                    Event new_event = Event(*bkpt1, nbp1, false);
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }
                if(getCircleEvent(nbp2, *bkpt2).y < directrix && check_clockwise(nbp2, *bkpt2)){
                    cout<<"Circle inserted using"<<endl;
                    auto x = *bkpt2;
                    nbp2.printit();
                    x.printit();
                    Event new_event = Event(nbp2, *bkpt2, false);
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }
                beachline.insert(breakpoint(bkpt1->f2, p0));
                beachline.insert(breakpoint(p0, bkpt2->f1));
            }
        }
    }

    void handle_circle_event(Event e){
        auto it = events.find(e);
        directrix = e.p.y;  
        if(it != events.end()){
            auto b1 = std::find(beachline.begin(), beachline.end(), it->b1);
            auto b2 = std::find(beachline.begin(), beachline.end(), it->b2);
            auto b0 = std::prev(b1);    
            auto b3 = std::next(b2); 
            
            auto b4 = breakpoint(b1->f1, b2->f2);
            auto a = *b1;
            auto k = *b2;
            std::cout<<"This is b1"<<endl;
            a.printit();
            std::cout<<"This is b2"<<endl;
            k.printit();
            cout<<"Circle event at "<<e.p.x<<" "<<e.p.y<<endl;

            if(!(b1 == beachline.begin())){
                if(b0->f2 == b1->f1 && b0->f1 == b1->f2){
                    cout<<"NO CIRCLE ABOVE"<<endl;
                }
                else{
                    cout<<"CIRCLE MAY BE ABOVE"<<endl;   
                    Event f = Event(*b0, *b1, false);
                    if(events.lower_bound(f) != events.end()){
                        cout<<"Circle found and erased ("<<f.p.x<<" "<<f.p.y<<")"<<endl;
                        events.erase(f);
                    }
                }
                if(getCircleEvent(*b0, b4).y <= directrix && check_clockwise(*b0, b4)){
                    cout<<"Circle inserted using"<<endl;
                    auto x = *b0;
                    x.printit();
                    b4.printit();
                    Event new_event = Event(*b0, b4, false);
                    cout<<"Inserting new circle event on the left"<<endl;
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }
            }

            if(!(b3 == beachline.end())){
                if(b2->f2 == b3->f1 && b2->f1 == b3->f2){
                    cout<<"NO CIRCLE ABOVE"<<endl;
                }
                else{
                    cout<<"CIRCLE MAY BE ABOVE"<<endl;   
                    Event f = Event(*b2, *b3, false);
                    if(events.lower_bound(f) != events.end()){
                        cout<<"Circle found and erased ("<<f.p.x<<" "<<f.p.y<<")"<<endl;
                        events.erase(f);
                        
                    }
                }
                if(getCircleEvent(b4, *b3).y <= directrix && check_clockwise(b4, *b3)){
                    cout<<"Circle inserted using"<<endl;
                    auto x = *b3;
                    b4.printit();
                    x.printit();
                    cout<<"Inserting new circle event on the right "<<endl;
                    Event new_event = Event(b4, *b3, false);
                    cout<<"("<<new_event.p.x<<" "<<new_event.p.y<<")"<<endl;
                    events.insert(new_event);
                }
                
            }
            beachline.erase(b1);
            beachline.erase(b2);
            beachline.insert(b4);
          
        }
        
    }
};


int main(){
    vector<Point> points = {Point(0, 10), Point(0.1, 5), Point(0.2, 3), Point(8, 0.5)};
    // vector<Point> points = {Point(1, 7), Point(5, 20), Point(50, 30), Point(7, -100)};
    Beachline b;
    for (auto p : points) {
        b.events.insert(Event(p, true));
    }
    std::cout<<b.events.size()<<endl;
    
    while(!b.events.empty()){
        Event e = *b.events.begin();
        directrix = e.p.y;
        if(e.is_site_event) {
            cout<<"----------"<<endl;
            cout<<"Site event ";
            cout<<"("<<e.p.x<<" "<<e.p.y<<")"<<endl;
            b.handle_site_event(e.p);
        }
        else{
            cout<<"Circle event";
            cout<<"("<<e.p.x<<" "<<e.p.y<<")"<<endl;
            b.handle_circle_event(e);
        }
        cout<<"Beachline"<<endl;
        for(auto i: b.beachline){i.printit();}
        std::cout<<endl;
        b.events.erase(b.events.begin());
        cout<<"event queue"<<endl;
        for(auto i: b.events){
            if(i.is_site_event) cout<<"Site event ";
            else cout<<"Circle event";
            cout<<"("<<i.p.x<<" "<<i.p.y<<")"<<endl;
        }
        cout<<"----------"<<endl;
    }
    
    cout<<endl;

    for(auto i: b.beachline){
        i.printit();}
    
}
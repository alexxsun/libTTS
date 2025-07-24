#ifndef MORSEINCIDENCEGRAPH
#define MORSEINCIDENCEGRAPH

#include "../iastar/simplicialcomplex.h"

class Node;

class Arc {
private:
    bool deleted;
    int label;

    Node* nBig;
    implicitS fromBig; // simplex of bigger dimension

    Node* nSmall;
    implicitS fromSmall; // simplex of sma

public:
    inline Arc(Node* u, implicitS firstn1, Node* d, implicitS firstn2)
        : nBig(u), fromBig(firstn1), nSmall(d), fromSmall(firstn2), label(1), deleted(false) {
    }

    inline Node* getNodeBig() { return nBig; }
    inline implicitS getFromBig() { return fromBig; }
    inline Node* getNodeSmall() { return nSmall; }
    inline implicitS getFromSmall() { return fromSmall; }

    inline Node* getOtherNode(Node* n) {
        if (n == nBig)
            return nSmall;
        else if (n == nSmall)
            return nBig;
        else
            return NULL;
    }

    inline int getLabel() { return label; }
    inline void setLabel(int l) { label = l; }
    inline void setDeleted() { deleted = true; }
    inline bool isDeleted() { return deleted; }
};

class Node {
private:
    implicitS csimpl;
    set<Arc*> arcsUp;
    set<Arc*> arcsDown;

public:
    inline Node(implicitS s) {
        csimpl = s;
        arcsUp = set<Arc*>();
        arcsDown = set<Arc*>();
    }

    inline implicitS& getCSimplex() { return csimpl; }
    inline set<Arc*>& getArcsUp() { return arcsUp; }
    inline set<Arc*>& getArcsDown() { return arcsDown; }

    inline void removeArc(Arc* a) {
        if (a->getNodeSmall()->getCSimplex() == csimpl)
            arcsUp.erase(a);
        else
            arcsDown.erase(a);
    }

    inline void addArc(Arc* a) {
        if (a->getNodeSmall()->getCSimplex() == csimpl)
            arcsUp.insert(a);
        else
            arcsDown.insert(a);
    }
};

class MIG {
    vector<set<Node*> > nodes;
    vector<set<Arc*> > arcs; // arcs connecting i-saddle and (i+1)-saddle are in position i

public:
    inline vector<set<Node*> >& getNodes() { return nodes; }
    inline vector<set<Arc*> >& getArcs() { return arcs; }

    inline MIG() {
    }

    inline MIG(int dim) {
        nodes = vector<set<Node*> >(dim, set<Node*>());
        arcs = vector<set<Arc*> >(dim, set<Arc*>());
    }

    inline Node* addNode(implicitS s) {
        Node* node = new Node(s);
        nodes[s.getDim()].insert(node);
        return node;
    }

    inline Arc* addArc(Node* n1, implicitS firstn1, Node* n2, implicitS firstn2) {
        Arc* a = new Arc(n1, firstn1, n2, firstn2);
        arcs[n2->getCSimplex().getDim()].insert(a);
        n1->addArc(a);
        n2->addArc(a);
        return a;
    }

    inline void removeArc(Arc* arc) {
        arcs[arc->getNodeSmall()->getCSimplex().getDim()].erase(arc);
        arc->getNodeBig()->removeArc(arc);
        arc->getNodeSmall()->removeArc(arc);
        arc->setDeleted();
    }

    inline void deleteNode(Node* node) {
        for (auto a : node->getArcsDown()) {
            arcs[a->getNodeSmall()->getCSimplex().getDim()].erase(a);
            a->getNodeSmall()->removeArc(a);
            a->setDeleted();
        }

        for (auto a : node->getArcsUp()) {
            arcs[a->getNodeSmall()->getCSimplex().getDim()].erase(a);
            a->getNodeBig()->removeArc(a);
            a->setDeleted();
        }

        nodes[node->getCSimplex().getDim()].erase(node);
        // delete node;
    }
};

#endif  // MORSEINCIDENCEGRAPH
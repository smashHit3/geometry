#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

using Int128 = __int128_t;

struct Node {
    long long x;
    long long y;
    long long value;
    long long sum;
    long long minimumX;
    long long maximumX;
    long long minimumY;
    long long maximumY;
    int size;
    int left;
    int right;
};

std::vector<Node> tree(1);
int root = 0;

void pull(int index) {
    Node& node = tree[index];
    node.sum = node.value;
    node.size = 1;
    node.minimumX = node.maximumX = node.x;
    node.minimumY = node.maximumY = node.y;
    for (const int child : {node.left, node.right}) {
        if (child == 0) {
            continue;
        }
        const Node& other = tree[child];
        node.sum += other.sum;
        node.size += other.size;
        node.minimumX = std::min(node.minimumX, other.minimumX);
        node.maximumX = std::max(node.maximumX, other.maximumX);
        node.minimumY = std::min(node.minimumY, other.minimumY);
        node.maximumY = std::max(node.maximumY, other.maximumY);
    }
}

int newNode(long long x, long long y, long long value) {
    tree.push_back({x, y, value, value, x, x, y, y, 1, 0, 0});
    return static_cast<int>(tree.size()) - 1;
}

void collect(int index, std::vector<int>& indices) {
    if (index == 0) {
        return;
    }
    indices.push_back(index);
    collect(tree[index].left, indices);
    collect(tree[index].right, indices);
}

int build(std::vector<int>& indices, int left, int right) {
    if (left >= right) {
        return 0;
    }

    long long minimumX = tree[indices[left]].x;
    long long maximumX = minimumX;
    long long minimumY = tree[indices[left]].y;
    long long maximumY = minimumY;
    for (int i = left + 1; i < right; ++i) {
        const Node& node = tree[indices[i]];
        minimumX = std::min(minimumX, node.x);
        maximumX = std::max(maximumX, node.x);
        minimumY = std::min(minimumY, node.y);
        maximumY = std::max(maximumY, node.y);
    }
    const int axis =
        static_cast<Int128>(maximumX) - minimumX >= static_cast<Int128>(maximumY) - minimumY
            ? 0
            : 1;
    const int middle = left + (right - left) / 2;
    std::nth_element(indices.begin() + left, indices.begin() + middle, indices.begin() + right,
                     [axis](int a, int b) {
                         const Node& first = tree[a];
                         const Node& second = tree[b];
                         const long long firstCoordinate = axis == 0 ? first.x : first.y;
                         const long long secondCoordinate = axis == 0 ? second.x : second.y;
                         if (firstCoordinate != secondCoordinate) {
                             return firstCoordinate < secondCoordinate;
                         }
                         const long long firstOther = axis == 0 ? first.y : first.x;
                         const long long secondOther = axis == 0 ? second.y : second.x;
                         return firstOther < secondOther;
                     });

    const int root = indices[middle];
    tree[root].left = build(indices, left, middle);
    tree[root].right = build(indices, middle + 1, right);
    pull(root);
    return root;
}

void rebuild(int& root) {
    std::vector<int> indices;
    indices.reserve(tree[root].size);
    collect(root, indices);
    root = build(indices, 0, static_cast<int>(indices.size()));
}

bool unbalanced(int index) {
    const Node& node = tree[index];
    const int largestChild = std::max(tree[node.left].size, tree[node.right].size);
    return largestChild * 4 > node.size * 3;
}

void insert(int& root, long long x, long long y, long long value) {
    if (root == 0) {
        root = newNode(x, y, value);
        return;
    }

    if (tree[root].x == x && tree[root].y == y) {
        tree[root].value += value;
        pull(root);
        return;
    }

    const Int128 width = static_cast<Int128>(tree[root].maximumX) - tree[root].minimumX;
    const Int128 height = static_cast<Int128>(tree[root].maximumY) - tree[root].minimumY;
    const bool splitX = width >= height;
    const bool goesLeft =
        splitX ? (x < tree[root].x || (x == tree[root].x && y < tree[root].y))
               : (y < tree[root].y || (y == tree[root].y && x < tree[root].x));
    if (goesLeft) {
        int child = tree[root].left;
        insert(child, x, y, value);
        tree[root].left = child;
    } else {
        int child = tree[root].right;
        insert(child, x, y, value);
        tree[root].right = child;
    }
    pull(root);
    if (unbalanced(root)) {
        rebuild(root);
    }
}

bool disjoint(const Node& node, long long x1, long long y1, long long x2, long long y2) {
    return node.maximumX < x1 || x2 < node.minimumX || node.maximumY < y1 ||
           y2 < node.minimumY;
}

bool contained(const Node& node, long long x1, long long y1, long long x2, long long y2) {
    return x1 <= node.minimumX && node.maximumX <= x2 && y1 <= node.minimumY &&
           node.maximumY <= y2;
}

long long rectangleSum(int root, long long x1, long long y1, long long x2, long long y2) {
    if (root == 0 || disjoint(tree[root], x1, y1, x2, y2)) {
        return 0;
    }
    const Node& node = tree[root];
    if (contained(node, x1, y1, x2, y2)) {
        return node.sum;
    }
    const long long own = x1 <= node.x && node.x <= x2 && y1 <= node.y && node.y <= y2
                              ? node.value
                              : 0;
    return own + rectangleSum(node.left, x1, y1, x2, y2) +
           rectangleSum(node.right, x1, y1, x2, y2);
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    tree.reserve(200000 + 1);
    long long lastAnswer = 0;
    int operation;
    while (std::cin >> operation && operation != 3) {
        if (operation == 0) {
            int boardSize;
            std::cin >> boardSize;
            continue;
        }
        if (operation == 1) {
            long long x;
            long long y;
            long long value;
            std::cin >> x >> y >> value;
            x ^= lastAnswer;
            y ^= lastAnswer;
            value ^= lastAnswer;
            insert(root, x, y, value);
        } else if (operation == 2) {
            long long x1;
            long long y1;
            long long x2;
            long long y2;
            std::cin >> x1 >> y1 >> x2 >> y2;
            x1 ^= lastAnswer;
            y1 ^= lastAnswer;
            x2 ^= lastAnswer;
            y2 ^= lastAnswer;
            if (x1 > x2) {
                std::swap(x1, x2);
            }
            if (y1 > y2) {
                std::swap(y1, y2);
            }
            lastAnswer = rectangleSum(root, x1, y1, x2, y2);
            std::cout << lastAnswer << '\n';
        }
    }
    return 0;
}

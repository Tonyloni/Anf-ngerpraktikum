#include <vector>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <iostream>

//T is the type of the element, K is the key, so the value that gets compared
template <typename T, typename K>
class Adr_PriorityQueue {
public:

    Adr_PriorityQueue() {}

    std::size_t size() const {
        return heap_.size();
    }

    bool isEmpty() const {
        return heap_.empty();
    }

    bool isElementOf(T element) const {
        return addressMap.find(element) != addressMap.end();
    }

    void insert(T element, K key) {
        Node newNode(element, key);
        heap_.push_back(newNode);
        std::size_t index = heap_.size() - 1;
        newNode.index = index;
        addressMap[element] = index;
        //std::cout << "addressMap[" << element << "]: " << addressMap[element] << std::endl;
        swiftUp(index);
    }

    T deleteMin() {
        if (heap_.empty()) {
            std::cout << "heap is empty" << std::endl;
        }
        T minElement = heap_[0].data;
        addressMap.erase(minElement);
        heap_[0].index = heap_.size() - 1;
        swapNodes(0, heap_.size()-1);
        heap_.pop_back();
        swiftDown(0);
        return minElement;
    }

    T top() {
        if (heap_.empty()) {
            std::cout << "heap is empty" << std::endl;
        }
        return heap_[0].data;
    }

    void updateKey(T element, K newKey) {


        std::size_t index = addressMap[element];
        //std::cout << "index: " << index << std::endl;
        Node& node = heap_[index];

        //std::cout << "newKey: " << newKey << std::endl;
        //std::cout << "node.key: " << node.key << std::endl;

        if (newKey < node.key) {
            node.key = newKey;
            swiftUp(index);
        } else if (newKey > node.key) {
           // std::cout << "HERE" << std::endl;
            node.key = newKey;
            swiftDown(index);
        }
    }

    K getKey(T element) {
        std::size_t index = addressMap[element];
        return heap_[index].key;
    }

private:
  struct Node {
    T data;
    K key;
    std::size_t index;

    Node(T data, K key)
      : data(data), key(key), index(0) {}

    bool operator<(const Node& other) const {
      return key < other.key;
    }
  };

  std::vector<Node> heap_;
  std::unordered_map<T, std::size_t> addressMap;

  void printMap() {
    for (int i = 1; i <= addressMap.size(); i++) {
        std::cout << "addressMap[" << i << "]: " << addressMap[i] << std::endl;
    }
    std::cout << std::endl;
  }

    void swapNodes(std::size_t index1, std::size_t index2) {
        addressMap[heap_[index1].data] = index2;
        addressMap[heap_[index2].data] = index1;
        std::swap(heap_[index1], heap_[index2]);
        //std::cout << "heap_[index1].data: " << heap_[index1].data << std::endl;
        //std::cout << "heap_[index2].data: " << heap_[index2].data << std::endl;
        
        //heap_[index1].index = index1;
        //heap_[index2].index = index2;
    }

  void swiftUp(std::size_t index) {
    while (index > 0) {
      std::size_t parentIndex = (index - 1) / 2;
      if (heap_[index] < heap_[parentIndex]) {
        swapNodes(index, parentIndex);
        index = parentIndex;
      } else {
        break;
      }
    }
  }

  void swiftDown(std::size_t index) {
    //printMap();
    std::size_t leftChild = 2 * index + 1;
    std::size_t rightChild = 2 * index + 2;
    std::size_t smallest = index;
    if (leftChild < heap_.size() && heap_[leftChild].key < heap_[smallest].key) {
      smallest = leftChild;
    }
    if (rightChild < heap_.size() && heap_[rightChild].key < heap_[smallest].key) {
      smallest = rightChild;
    }
    if (smallest != index) {
        //std::cout << "Index: " << index << ", smallest: " << smallest << std::endl;
      swapNodes(index, smallest);
      swiftDown(smallest);
    }
  }
};


#include <iostream>
#include <unordered_set>
#include <vector>
#include <random>
#include <cassert>

const std::unordered_set<std::string> unary_operators = {"~", "log", "ln", "exp", "cos", "sin", "sqrt", "asin", "arcsin", "acos", "arccos", "tanh", "sech", "abs"};
const std::unordered_set<std::string> binary_operators =  {"+", "-", "*", "/", "^"};
const std::unordered_set<std::string> leaf_nodes = {"x", "y", "z", "0", "1", "2", "4"};
std::vector<std::string> unary_ops(unary_operators.begin(), unary_operators.end()), binary_ops(binary_operators.begin(), binary_operators.end()), leaves(leaf_nodes.begin(), leaf_nodes.end());
int N = 5;

bool is_unary(const std::string& token)
{
    return (unary_operators.find(token) != unary_operators.end());
}

bool is_binary(const std::string& token)
{
    return (binary_operators.find(token) != binary_operators.end());
}

bool is_leaf(const std::string& token)
{
    return (!is_unary(token) && !is_binary(token));
}

std::pair<int, bool> getRPNdepth(const std::vector<std::string>& expression)
{
    if (expression.empty())
    {
        return std::make_pair(0, false);
    }
    
    std::vector<int> stack;
    bool complete = true;
    for (size_t i = 0; i < expression.size(); i++)
    {
        if (is_unary(expression[i]))
        {
            stack.back() += 1;
        }
        else if (is_binary(expression[i]))
        {
            int op2 = stack.back();
            stack.pop_back();
            int op1 = stack.back();
            stack.pop_back();
            stack.push_back(std::max(op1, op2) + 1);
        }
        else //leaf
        {
            stack.push_back(1);
        }
    }

    while (stack.size() > 1)
    {
        int op2 = stack.back();
        stack.pop_back();
        int op1 = stack.back();
        stack.pop_back();
        stack.push_back(std::max(op1, op2) + 1);
        complete = false;
    }

    /*
     e.g., assume stack = {1, 2, 3, 4, 5}, then:
     {1, 2, 3, 4, 5}
     {1, 2, 3, 6}
     {1, 2, 7}
     {1, 8}
     {9}
     */

    return std::make_pair(stack.back() - 1, complete);
}

int numBinaryOps(const std::vector<std::string>& expression)
{
    int count = 0;
    for (const auto& i: expression)
    {
        if (is_binary(i))
        {
            count++;
        }
    }
    return count;
}

int numLeaves(const std::vector<std::string>& expression)
{
    int count = 0;
    for (const auto& i: expression)
    {
        if (is_leaf(i))
        {
            count++;
        }
    }
    return count;
}

std::vector<std::string> get_legal_moves(std::vector<std::string>& expression)
{
    if (expression.empty()) //At the beginning, expression is empty, so the only legal moves are the leaves
    {
        return leaves;
    }
    int num_binary = numBinaryOps(expression);
    int num_leaves = numLeaves(expression);

    bool una_allowed = false, bin_allowed = (num_binary != num_leaves - 1), leaf_allowed = false;
    if (unary_operators.size() > 0)
    {
        expression.push_back(unary_ops[0]);
        una_allowed = ((num_leaves >= 1) && (getRPNdepth(expression).first <= N));
    }

    expression[expression.size() - 1] = leaves[0];
    leaf_allowed = (getRPNdepth(expression).first <= N);

    expression.pop_back();

    std::vector<std::string> allowed_tokens;
    if (una_allowed)
    {
        for (const auto& i: unary_ops)
        {
            allowed_tokens.push_back(i);
        }
    }
    if (bin_allowed)
    {
        for (const auto& i: binary_ops)
        {
            allowed_tokens.push_back(i);
        }
    }
    if (leaf_allowed)
    {
        for (const auto& i: leaf_nodes)
        {
            allowed_tokens.push_back(i);
        }
    }
    return allowed_tokens;
}

bool complete(const std::vector<std::string>& expression)
{
    auto [depth, complete] = getRPNdepth(expression);
    if (complete && (depth == N))
    {
        return true;
    }
    return false;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    for (const T& i: vec)
    {
        os << i << ' ';
    }
    return os;
}

std::vector<std::string> generate_func()
{
    std::vector<std::string> expr;
    std::vector<std::string> legal_moves;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    while (!complete(expr))
    {
        legal_moves = get_legal_moves(expr);
        std::uniform_int_distribution<int> distribution(0, legal_moves.size() - 1);
        expr.push_back(legal_moves[distribution(generator)]);
    }
    return expr;
}

std::string to_infix(const std::vector<std::string>& expression)
{
    std::stack<std::string> stack;
    std::string result, token;
    int sz = expression.size();
    for (int i = 0; i < sz; i++)
    {
        token = expression[i];
        if (is_leaf(token)) // leaf
        {
            stack.push(token);
        }
        else if (is_unary(token)) // Unary operator
        {
            std::string operand = stack.top();
            stack.pop();
            result = token + "(" + operand + ")";
            stack.push(result);
        }
        else // binary operator
        {
            std::string right_operand = stack.top();
            stack.pop();
            std::string left_operand = stack.top();
            stack.pop();
            result = "(" + left_operand + " " + token + " " + right_operand + ")";
            stack.push(result);
        }
    }

    return stack.top();
}

std::string to_infix_latex(const std::vector<std::string>& expression)
{
    std::stack<std::string> stack;
    std::string result, token;
    int sz = expression.size();
    for (int i = 0; i < sz; i++)
    {
        token = expression[i];
        if (is_leaf(token)) // leaf
        {
            stack.push(token);
        }
        else if (is_unary(token)) // Unary operator
        {
            std::string operand = stack.top();
            
            stack.pop();
            result = "\\mathrm{" + token + "}\\left(" + operand + "\\right)";
            stack.push(result);
        }
        else // binary operator
        {
            std::string right_operand = stack.top();
            stack.pop();
            std::string left_operand = stack.top();
            stack.pop();
            result = "\\left(" + left_operand + " " + token + " " + right_operand + "\\right)";
            stack.push(result);
        }
    }

    return stack.top();
}

int main(int argc, char** argv)
{
    assert(argc == 2);
    int num = std::stoi(argv[1]);
    for (int i = 0; i < num; i++)
    {
        std::cout << to_infix_latex(generate_func()) << '\n';
    }
}

//g++ -std=c++20 -o fixedExprDepthGenerator fixedExprDepthGenerator.cpp -L/opt/homebrew/Cellar/boost/1.84.0 -I/opt/homebrew/Cellar/boost/1.84.0/include
